#!/usr/bin/env python3
"""
compare_trees.py — compare legacy (scalar-per-row) vs SoA (vector-per-event)
TriggerAnaTree output trees.

Usage:
    python compare_trees.py <rootfile> [--rtol 1e-6] [--atol 1e-9] [--verbose]
"""

import argparse
import sys
import numpy as np
import uproot
from rich.console import Console
from rich.table import Table
from rich import print as rprint

console = Console()


# ---------------------------------------------------------------------------
# Tree pair definitions
# ---------------------------------------------------------------------------

# Fixed top-level pairs: (legacy_path, soa_path)
FIXED_PAIRS = [
    ("mctruths",     "mctruths_2g"),
    ("mcneutrinos",  "mcneutrinos_2g"),
    ("mcparticles",  "mcparticles_2g"),
    ("simides",      "simides_2g"),
]

# Scalar-per-event pairs (both trees have one row per event)
SCALAR_PAIRS = [
    ("event_summary", "event_summary_2g"),
]

# Subdirectories to scan for auto-discovered TP/TA/TC pairs
SCAN_DIRS = ["TriggerPrimitives"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def event_key(run, subrun, event):
    return (int(run), int(subrun), int(event))


def is_float_branch(arr):
    return np.issubdtype(arr.dtype, np.floating)


def sort_index(data: dict, sort_branch: str):
    """Return argsort on the first available sort_branch."""
    for b in sort_branch if isinstance(sort_branch, list) else [sort_branch]:
        if b in data and len(data[b]):
            return np.argsort(data[b], kind="stable")
    # Fallback: identity
    n = next(iter(data.values()))
    return np.arange(len(n))


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------

def load_legacy_tree(tree):
    """
    Load a legacy tree (one row per object).
    Returns dict: event_key -> {branch: np.ndarray}
    """
    branches = [b for b in tree.keys() if b not in ("event", "run", "subrun")]
    meta = tree.arrays(["run", "subrun", "event"], library="np")
    data = tree.arrays(branches, library="np")

    grouped = {}
    for i in range(len(meta["event"])):
        key = event_key(meta["run"][i], meta["subrun"][i], meta["event"][i])
        if key not in grouped:
            grouped[key] = {b: [] for b in branches}
        for b in branches:
            grouped[key][b].append(data[b][i])

    # Convert lists to numpy arrays
    result = {}
    for key, bd in grouped.items():
        result[key] = {b: np.array(v) for b, v in bd.items()}
    return result


def load_soa_tree(tree):
    """
    Load a SoA tree (one row per event, fields are vectors).
    Returns dict: event_key -> {branch: np.ndarray}
    """
    branches = [b for b in tree.keys() if b not in ("event", "run", "subrun")]
    meta = tree.arrays(["run", "subrun", "event"], library="np")
    # SoA branches are std::vector<T> — use awkward or iterate
    data = tree.arrays(branches, library="ak")

    result = {}
    n_events = len(meta["event"])
    for i in range(n_events):
        key = event_key(meta["run"][i], meta["subrun"][i], meta["event"][i])
        result[key] = {}
        for b in branches:
            val = data[b][i]
            # Convert awkward array to numpy
            try:
                result[key][b] = np.asarray(val)
            except Exception:
                result[key][b] = np.array(val.tolist())
    return result


def load_scalar_tree(tree):
    """
    Load a scalar-per-event tree (event_summary style).
    Returns dict: event_key -> {branch: scalar}
    """
    branches = [b for b in tree.keys() if b not in ("event", "run", "subrun")]
    meta = tree.arrays(["run", "subrun", "event"], library="np")
    data = tree.arrays(branches, library="np")

    result = {}
    for i in range(len(meta["event"])):
        key = event_key(meta["run"][i], meta["subrun"][i], meta["event"][i])
        result[key] = {b: data[b][i] for b in branches}
    return result


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

# Preferred sort branch per tree type (tried in order)
SORT_BRANCHES = {
    "mctruths":    ["block_id", "truth_track_id", "pdg"],
    "mcneutrinos": ["block_id", "nupdg"],
    "mcparticles": ["g4_track_id", "pdg"],
    "simides":     ["channel", "timestamp", "trackID"],
}
DEFAULT_SORT = ["channel", "block_id", "g4_track_id", "trackID", "time_start"]


def compare_object_trees(name_legacy, tree_legacy, tree_soa, rtol, atol, verbose):
    """Compare a legacy (scalar-per-row) tree vs a SoA (vector-per-event) tree."""
    print(f"\n--- Comparing: {name_legacy}  vs  {name_legacy}_2g ---")

    leg = load_legacy_tree(tree_legacy)
    soa = load_soa_tree(tree_soa)

    all_keys = sorted(set(leg) | set(soa))
    if set(leg) != set(soa):
        only_leg = sorted(set(leg) - set(soa))
        only_soa = sorted(set(soa) - set(leg))
        if only_leg:
            print(f"  [WARN] Events only in legacy: {only_leg[:5]}")
        if only_soa:
            print(f"  [WARN] Events only in SoA:    {only_soa[:5]}")

    # Determine sort branches for this tree
    sort_branches = SORT_BRANCHES.get(name_legacy.rstrip("/").split("/")[-1], DEFAULT_SORT)

    n_mismatch = 0
    branch_mismatches = {}

    for key in all_keys:
        if key not in leg or key not in soa:
            continue
        ld = leg[key]
        sd = soa[key]

        # Check object count
        first_b = next(iter(ld))
        n_leg = len(ld[first_b])
        n_soa = len(sd.get(first_b, []))
        if n_leg != n_soa:
            n_mismatch += 1
            if verbose:
                print(f"  Event {key}: object count mismatch: legacy={n_leg}, soa={n_soa}")
            continue

        # Sort both by stable key
        idx_leg = sort_index(ld, sort_branches)
        idx_soa = sort_index(sd, sort_branches)

        event_ok = True
        all_sorted = {}   # branch -> (lv_sorted, sv_sorted)
        diff_masks = {}   # branch -> diff_mask  (only failing branches)

        for b in ld:
            if b not in sd:
                continue
            lv = ld[b][idx_leg]
            sv = sd[b][idx_soa]
            all_sorted[b] = (lv, sv)

            if lv.dtype.kind in ('U', 'S', 'O'):  # strings
                ok = np.array_equal(lv, sv)
                diff_mask = lv != sv
            elif is_float_branch(lv):
                ok = np.allclose(lv, sv, rtol=rtol, atol=atol, equal_nan=True)
                diff_mask = ~np.isclose(lv, sv, rtol=rtol, atol=atol, equal_nan=True)
            else:
                ok = np.array_equal(lv, sv)
                diff_mask = lv != sv

            if not ok:
                event_ok = False
                branch_mismatches[b] = branch_mismatches.get(b, 0) + 1
                diff_masks[b] = diff_mask

        if not event_ok:
            n_mismatch += 1
            if verbose and diff_masks:
                # Union of all differing row indexes across branches
                diff_mask_all = np.zeros(n_leg, dtype=bool)
                for mask in diff_masks.values():
                    diff_mask_all |= mask
                row_idx = np.where(diff_mask_all)[0]

                table = Table(
                    title=f"Event {key} — {len(row_idx)} mismatching object(s) of {n_leg}",
                    show_lines=False,
                    highlight=True,
                )
                table.add_column("idx", justify="right", style="dim")
                for b in sorted(all_sorted):
                    bad_col = b in diff_masks
                    style = "bold red" if bad_col else ""
                    table.add_column(f"{b} [L]", justify="right", style=style)
                    table.add_column(f"{b} [S]", justify="right", style=style)

                for i in row_idx:
                    cells = [str(i)]
                    for b in sorted(all_sorted):
                        lv, sv = all_sorted[b]
                        bad_cell = b in diff_masks and diff_masks[b][i]
                        lstr = f"[bold red]{lv[i]}*[/bold red]" if bad_cell else str(lv[i])
                        sstr = str(sv[i])
                        cells += [lstr, sstr]
                    table.add_row(*cells)

                console.print(table)

    status = "PASS" if n_mismatch == 0 else "FAIL"
    print(f"  Result: {status}  ({n_mismatch}/{len(all_keys)} events with mismatches)")
    if branch_mismatches:
        for b, cnt in sorted(branch_mismatches.items(), key=lambda x: -x[1]):
            print(f"    branch '{b}': {cnt} mismatching events")
    return n_mismatch == 0


def compare_scalar_trees(name, tree_legacy, tree_soa, rtol, atol, verbose):
    """Compare two scalar-per-event trees (event_summary style)."""
    print(f"\n--- Comparing: {name}  vs  {name}_2g ---")

    leg = load_scalar_tree(tree_legacy)
    soa = load_scalar_tree(tree_soa)

    all_keys = sorted(set(leg) | set(soa))
    n_mismatch = 0

    for key in all_keys:
        if key not in leg or key not in soa:
            n_mismatch += 1
            continue
        ld = leg[key]
        sd = soa[key]
        event_ok = True
        scalar_diffs = {}
        for b in ld:
            if b not in sd:
                continue
            lv, sv = ld[b], sd[b]
            if isinstance(lv, float):
                ok = np.isclose(lv, sv, rtol=rtol, atol=atol)
            else:
                ok = (lv == sv)
            if not ok:
                event_ok = False
                scalar_diffs[b] = (lv, sv)
        if not event_ok:
            n_mismatch += 1
            if verbose and scalar_diffs:
                table = Table(title=f"Event {key}", show_lines=False, highlight=True)
                table.add_column("branch", style="dim")
                table.add_column("legacy [L]", justify="right", style="bold red")
                table.add_column("soa [S]",    justify="right")
                for b, (lv, sv) in sorted(scalar_diffs.items()):
                    table.add_row(b, str(lv), str(sv))
                console.print(table)

    status = "PASS" if n_mismatch == 0 else "FAIL"
    print(f"  Result: {status}  ({n_mismatch}/{len(all_keys)} events with mismatches)")
    return n_mismatch == 0


# ---------------------------------------------------------------------------
# Auto-discover TP/TA/TC tree pairs in subdirectories
# ---------------------------------------------------------------------------

def discover_pairs(f, subdir):
    """
    Find all trees in f[subdir] whose name does NOT end in '_2g',
    then check for a sibling ending in '_2g'.
    Returns list of (legacy_path, soa_path).
    """
    pairs = []
    try:
        keys = [k for k in f[subdir].keys(cycle=False)]
    except Exception:
        return pairs

    names = set()
    for k in keys:
        # uproot keys may include ';1' cycle suffix
        name = k.split(";")[0]
        names.add(name)

    for name in sorted(names):
        if name.endswith("_2g"):
            continue
        soa_name = name + "_2g"
        if soa_name in names:
            pairs.append((f"{subdir}/{name}", f"{subdir}/{soa_name}"))
    return pairs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compare legacy vs SoA TriggerAnaTree output.")
    parser.add_argument("rootfile", help="Path to the ROOT file to compare")
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance for float comparison")
    parser.add_argument("--atol", type=float, default=1e-9, help="Absolute tolerance for float comparison")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print per-branch mismatch details")
    args = parser.parse_args()

    print(f"Opening: {args.rootfile}")
    f = uproot.open(args.rootfile)

    all_pass = True

    # Fixed object-level pairs
    for leg_name, soa_name in FIXED_PAIRS:
        leg_name = 'triggerAna/' + leg_name
        soa_name = 'triggerAna/' + soa_name
        if leg_name not in f or soa_name not in f:
            print(f"\n[SKIP] {leg_name} or {soa_name} not found in file")
            continue
        ok = compare_object_trees(
            leg_name, f[leg_name], f[soa_name],
            args.rtol, args.atol, args.verbose
        )
        all_pass = all_pass and ok

    # Scalar-per-event pairs
    for leg_name, soa_name in SCALAR_PAIRS:
        leg_name = 'triggerAna/' + leg_name
        soa_name = 'triggerAna/' + soa_name
        if leg_name not in f or soa_name not in f:
            print(f"\n[SKIP] {leg_name} or {soa_name} not found in file")
            continue
        ok = compare_scalar_trees(
            leg_name, f[leg_name], f[soa_name],
            args.rtol, args.atol, args.verbose
        )
        all_pass = all_pass and ok

    # Auto-discovered pairs under subdirectories
    for subdir in SCAN_DIRS:
        subdir = 'triggerAna/' + subdir

        pairs = discover_pairs(f, subdir)
        if not pairs:
            print(f"\n[INFO] No tree pairs found under '{subdir}'")
        for leg_path, soa_path in pairs:
            ok = compare_object_trees(
                leg_path, f[leg_path], f[soa_path],
                args.rtol, args.atol, args.verbose
            )
            all_pass = all_pass and ok

    print("\n" + "=" * 50)
    print(f"Overall result: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 50)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()

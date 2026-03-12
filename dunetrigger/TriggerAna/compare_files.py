#!/usr/bin/env python3
"""
compare_files.py — compare TriggerAnaTree output across two ROOT files.

Both files are expected to use the SoA (vector-per-event) format with
identical tree names.

Usage:
    python compare_files.py <file_a> <file_b> [--rtol 1e-6] [--atol 1e-9] [--verbose]
"""

import argparse
import sys
import numpy as np
import uproot
from rich.console import Console
from rich.table import Table

console = Console()


# ---------------------------------------------------------------------------
# Tree definitions
# ---------------------------------------------------------------------------

# Top-level SoA object trees (one row per event, vector branches)
SOA_TREES = [
    "mctruths",
    "mcneutrinos",
    "mcparticles",
    "simides",
]

# Scalar-per-event trees
SCALAR_TREES = [
    "event_summary",
]

# Subdirectories to scan for auto-discovered TP trees
SCAN_DIRS = ["TriggerPrimitives"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def event_key(run, subrun, event):
    return (int(run), int(subrun), int(event))


def is_float_branch(arr):
    return np.issubdtype(arr.dtype, np.floating)


def sort_index(data: dict, sort_branches):
    """Return argsort on the first available sort branch."""
    for b in sort_branches if isinstance(sort_branches, list) else [sort_branches]:
        if b in data and len(data[b]):
            return np.argsort(data[b], kind="stable")
    n = len(next(iter(data.values())))
    return np.arange(n)


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------

def load_soa_tree(tree):
    """
    Load a tree, auto-detecting format:
      - SoA (one row per event, vector branches): each row maps directly to one event key.
      - Legacy (one row per object, scalar branches): rows are grouped by event key.
    Returns dict: event_key -> {branch: np.ndarray}
    """
    branches = [b for b in tree.keys() if b not in ("event", "run", "subrun")]
    meta = tree.arrays(["run", "subrun", "event"], library="np")
    if not branches or not len(meta["event"]):
        return {}
    data = tree.arrays(branches, library="ak")

    # Detect format: if the first branch of the first row is a scalar, it's a legacy tree.
    sample = np.asarray(data[branches[0]][0])
    is_legacy = sample.ndim == 0

    if is_legacy:
        # Accumulate scalar rows per event key (same as load_legacy_tree)
        grouped = {}
        for i in range(len(meta["event"])):
            key = event_key(meta["run"][i], meta["subrun"][i], meta["event"][i])
            if key not in grouped:
                grouped[key] = {b: [] for b in branches}
            for b in branches:
                grouped[key][b].append(np.asarray(data[b][i]).item())
        return {key: {b: np.array(v) for b, v in bd.items()} for key, bd in grouped.items()}
    else:
        # SoA: one row per event, each branch value is already an array
        result = {}
        for i in range(len(meta["event"])):
            key = event_key(meta["run"][i], meta["subrun"][i], meta["event"][i])
            result[key] = {}
            for b in branches:
                val = data[b][i]
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

SORT_BRANCHES = {
    "mctruths":    ["block_id", "truth_track_id", "pdg"],
    "mcneutrinos": ["block_id", "nupdg"],
    "mcparticles": ["g4_track_id", "pdg"],
    "simides":     ["channel", "timestamp", "trackID"],
}
DEFAULT_SORT = ["channel", "block_id", "g4_track_id", "trackID", "time_start"]


def compare_soa_trees(name, tree_a, tree_b, rtol, atol, verbose):
    """Compare two SoA trees (vector-per-event) with the same name from different files."""
    print(f"\n--- Comparing: {name} ---")

    da = load_soa_tree(tree_a)
    db = load_soa_tree(tree_b)

    all_keys = sorted(set(da) | set(db))
    if set(da) != set(db):
        only_a = sorted(set(da) - set(db))
        only_b = sorted(set(db) - set(da))
        if only_a:
            print(f"  [WARN] Events only in file A: {only_a[:5]}")
        if only_b:
            print(f"  [WARN] Events only in file B: {only_b[:5]}")

    sort_branches = SORT_BRANCHES.get(name.rstrip("/").split("/")[-1], DEFAULT_SORT)

    n_mismatch = 0
    branch_mismatches = {}

    for key in all_keys:
        if key not in da or key not in db:
            n_mismatch += 1
            continue

        ad = da[key]
        bd = db[key]

        common_branches = sorted(set(ad) & set(bd))
        if not common_branches:
            if verbose:
                print(f"  Event {key}: no common branches between files")
            n_mismatch += 1
            continue

        first_b = common_branches[0]
        n_a = len(ad[first_b])
        n_b = len(bd[first_b])
        if n_a != n_b:
            n_mismatch += 1
            if verbose:
                print(f"  Event {key}: object count mismatch: A={n_a}, B={n_b}")
            continue

        idx_a = sort_index(ad, sort_branches)
        idx_b = sort_index(bd, sort_branches)

        event_ok = True
        all_sorted = {}
        diff_masks = {}

        for b in common_branches:
            av = ad[b][idx_a]
            bv = bd[b][idx_b]
            all_sorted[b] = (av, bv)

            if av.dtype.kind in ('U', 'S', 'O'):
                ok = np.array_equal(av, bv)
                diff_mask = av != bv
            elif is_float_branch(av):
                ok = np.allclose(av, bv, rtol=rtol, atol=atol, equal_nan=True)
                diff_mask = ~np.isclose(av, bv, rtol=rtol, atol=atol, equal_nan=True)
            else:
                ok = np.array_equal(av, bv)
                diff_mask = av != bv

            if not ok:
                event_ok = False
                branch_mismatches[b] = branch_mismatches.get(b, 0) + 1
                diff_masks[b] = diff_mask

        if not event_ok:
            n_mismatch += 1
            if verbose and diff_masks:
                diff_mask_all = np.zeros(n_a, dtype=bool)
                for mask in diff_masks.values():
                    diff_mask_all |= mask
                row_idx = np.where(diff_mask_all)[0]

                table = Table(
                    title=f"Event {key} — {len(row_idx)} mismatching object(s) of {n_a}",
                    show_lines=False,
                    highlight=True,
                )
                table.add_column("idx", justify="right", style="dim")
                for b in sorted(all_sorted):
                    bad_col = b in diff_masks
                    style = "bold red" if bad_col else ""
                    table.add_column(f"{b} [A]", justify="right", style=style)
                    table.add_column(f"{b} [B]", justify="right", style=style)

                for i in row_idx:
                    cells = [str(i)]
                    for b in sorted(all_sorted):
                        av, bv = all_sorted[b]
                        bad_cell = b in diff_masks and diff_masks[b][i]
                        astr = f"[bold red]{av[i]}*[/bold red]" if bad_cell else str(av[i])
                        cells += [astr, str(bv[i])]
                    table.add_row(*cells)

                console.print(table)

    status = "PASS" if n_mismatch == 0 else "FAIL"
    print(f"  Result: {status}  ({n_mismatch}/{len(all_keys)} events with mismatches)")
    if branch_mismatches:
        for b, cnt in sorted(branch_mismatches.items(), key=lambda x: -x[1]):
            print(f"    branch '{b}': {cnt} mismatching events")
    return n_mismatch == 0


def compare_scalar_trees(name, tree_a, tree_b, rtol, atol, verbose):
    """Compare two scalar-per-event trees from different files."""
    print(f"\n--- Comparing: {name} ---")

    da = load_scalar_tree(tree_a)
    db = load_scalar_tree(tree_b)

    all_keys = sorted(set(da) | set(db))
    n_mismatch = 0

    for key in all_keys:
        if key not in da or key not in db:
            n_mismatch += 1
            continue
        ad = da[key]
        bd = db[key]
        event_ok = True
        scalar_diffs = {}
        for b in ad:
            if b not in bd:
                continue
            av, bv = ad[b], bd[b]
            ok = np.isclose(av, bv, rtol=rtol, atol=atol) if isinstance(av, float) else (av == bv)
            if not ok:
                event_ok = False
                scalar_diffs[b] = (av, bv)
        if not event_ok:
            n_mismatch += 1
            if verbose and scalar_diffs:
                table = Table(title=f"Event {key}", show_lines=False, highlight=True)
                table.add_column("branch", style="dim")
                table.add_column("file A", justify="right", style="bold red")
                table.add_column("file B", justify="right")
                for b, (av, bv) in sorted(scalar_diffs.items()):
                    table.add_row(b, str(av), str(bv))
                console.print(table)

    status = "PASS" if n_mismatch == 0 else "FAIL"
    print(f"  Result: {status}  ({n_mismatch}/{len(all_keys)} events with mismatches)")
    return n_mismatch == 0


# ---------------------------------------------------------------------------
# Auto-discover trees in subdirectories
# ---------------------------------------------------------------------------

def discover_trees(f, subdir):
    """
    Find all tree names in f[subdir] (stripping uproot cycle suffixes).
    Returns list of full paths.
    """
    paths = []
    try:
        keys = list(f[subdir].keys(cycle=False))
    except Exception:
        return paths
    for k in keys:
        paths.append(f"{subdir}/{k.split(';')[0]}")
    return sorted(set(paths))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compare TriggerAnaTree output across two ROOT files.")
    parser.add_argument("file_a", help="First ROOT file")
    parser.add_argument("file_b", help="Second ROOT file")
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance for float comparison")
    parser.add_argument("--atol", type=float, default=1e-9, help="Absolute tolerance for float comparison")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print per-branch mismatch details")
    args = parser.parse_args()

    print(f"File A: {args.file_a}")
    print(f"File B: {args.file_b}")
    fa = uproot.open(args.file_a)
    fb = uproot.open(args.file_b)

    all_pass = True

    # Fixed SoA object trees
    for name in SOA_TREES:
        path = f"triggerAna/{name}"
        if path not in fa or path not in fb:
            print(f"\n[SKIP] {path} not found in one or both files")
            continue
        ok = compare_soa_trees(path, fa[path], fb[path], args.rtol, args.atol, args.verbose)
        all_pass = all_pass and ok

    # Scalar-per-event trees
    for name in SCALAR_TREES:
        path = f"triggerAna/{name}"
        if path not in fa or path not in fb:
            print(f"\n[SKIP] {path} not found in one or both files")
            continue
        ok = compare_scalar_trees(path, fa[path], fb[path], args.rtol, args.atol, args.verbose)
        all_pass = all_pass and ok

    # Auto-discovered trees under subdirectories
    for subdir in SCAN_DIRS:
        subdir_path = f"triggerAna/{subdir}"
        trees_a = set(discover_trees(fa, subdir_path))
        trees_b = set(discover_trees(fb, subdir_path))
        common = sorted(trees_a & trees_b)
        only_a = sorted(trees_a - trees_b)
        only_b = sorted(trees_b - trees_a)

        if only_a:
            print(f"\n[WARN] Trees only in file A under {subdir_path}: {only_a}")
        if only_b:
            print(f"\n[WARN] Trees only in file B under {subdir_path}: {only_b}")
        if not common:
            print(f"\n[INFO] No common trees found under '{subdir_path}'")

        for path in common:
            ok = compare_soa_trees(path, fa[path], fb[path], args.rtol, args.atol, args.verbose)
            all_pass = all_pass and ok

    print("\n" + "=" * 50)
    print(f"Overall result: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 50)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()

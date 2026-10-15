#!/usr/bin/env python
"""Compare all TTrees (and the `info` TNamed JSON) in two TriggerAnaTree outputs.

usage: python compare_trees.py baseline.root new.root

PyROOT port of the uproot script in TriggerAnaTree_refactor_plan.md.
Objects are fetched by name without a cycle, so ROOT returns the latest cycle.
Exit code: 0 if identical, 1 otherwise.
"""
import json
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)


def walk(d, prefix=""):
    """Return {path: TTree} and {path: TNamed} for d and its subdirectories."""
    trees, nameds = {}, {}
    names = sorted({k.GetName() for k in d.GetListOfKeys()})  # one per name, any cycle
    for name in names:
        obj = d.Get(name)  # no cycle -> latest
        full = prefix + name
        if obj.InheritsFrom("TTree"):
            trees[full] = obj
        elif obj.InheritsFrom("TDirectory"):
            t, n = walk(obj, full + "/")
            trees.update(t)
            nameds.update(n)
        elif obj.InheritsFrom("TNamed") and obj.ClassName() == "TNamed":
            nameds[full] = obj
    return trees, nameds


def branch_names(tree):
    return {b.GetName() for b in tree.GetListOfBranches()}


def normalise(x):
    """Turn one entry (scalar, RVec, std::string, vector<string>) into a comparable value."""
    try:
        a = np.asarray(x)
        if a.dtype != object:
            return a
    except Exception:
        pass
    if isinstance(x, (str, bytes)):
        return str(x)
    try:
        return [str(s) for s in x]
    except TypeError:
        return str(x)


def same_column(xa, xb):
    if len(xa) != len(xb):
        return False
    if xa.dtype != object and xb.dtype != object:
        if xa.dtype.kind == "f" or xb.dtype.kind == "f":
            return np.array_equal(xa, xb, equal_nan=True)
        return np.array_equal(xa, xb)
    for u, v in zip(xa, xb):
        nu, nv = normalise(u), normalise(v)
        if isinstance(nu, np.ndarray) and isinstance(nv, np.ndarray):
            eq_nan = nu.dtype.kind == "f" and nv.dtype.kind == "f"
            if not np.array_equal(nu, nv, equal_nan=eq_nan):
                return False
        elif isinstance(nu, np.ndarray) or isinstance(nv, np.ndarray):
            return False
        elif nu != nv:
            return False
    return True


def read_columns(tree, columns):
    """Read columns with one event loop; fall back to one loop per column on failure."""
    if not columns:
        return {}, []
    try:
        return ROOT.RDataFrame(tree).AsNumpy(columns), []
    except Exception:
        data, bad = {}, []
        for c in columns:
            try:
                data.update(ROOT.RDataFrame(tree).AsNumpy([c]))
            except Exception:
                bad.append(c)
        return data, bad


def main(path_a, path_b):
    fa, fb = ROOT.TFile.Open(path_a), ROOT.TFile.Open(path_b)
    if not fa or fa.IsZombie() or not fb or fb.IsZombie():
        print("ERROR: cannot open input files")
        return 2
    print("OLD:", path_a)
    print("NEW:", path_b)
    ta, na = walk(fa)
    tb, nb = walk(fb)
    print("found", len(ta), "/", len(tb), "trees and", len(na), "/", len(nb), "TNamed (old / new)")
    diff = 0

    tree_names = sorted(set(ta) | set(tb))
    for i, name in enumerate(tree_names, 1):
        print("[%d/%d] %s" % (i, len(tree_names), name), flush=True)
        if name not in tb:
            print("ONLY IN OLD :", name); diff += 1; continue
        if name not in ta:
            print("ONLY IN NEW :", name); diff += 1; continue
        a, b = ta[name], tb[name]
        diff_before = diff
        if a.GetEntries() != b.GetEntries():
            print("ENTRIES     :", name, a.GetEntries(), "->", b.GetEntries()); diff += 1
        ka, kb = branch_names(a), branch_names(b)
        for br in sorted(ka - kb):
            print("BRANCH ONLY IN OLD:", name, br); diff += 1
        for br in sorted(kb - ka):
            print("BRANCH ONLY IN NEW:", name, br); diff += 1
        common = sorted(ka & kb)
        print("      reading", len(common), "common branches,", a.GetEntries(), "entries", flush=True)
        da, bad_a = read_columns(a, common)
        db, bad_b = read_columns(b, common)
        for br in common:
            if br in bad_a or br in bad_b:
                print("UNREADABLE  :", name, br); diff += 1
            elif not same_column(da[br], db[br]):
                print("CONTENT     :", name, br); diff += 1
        print("      ok" if diff == diff_before else "      %d difference(s)" % (diff - diff_before))

    print("comparing TNamed objects", flush=True)
    for name in sorted(set(na) | set(nb)):
        if name not in nb:
            print("TNAMED ONLY IN OLD:", name); diff += 1; continue
        if name not in na:
            print("TNAMED ONLY IN NEW:", name); diff += 1; continue
        sa, sb = na[name].GetTitle(), nb[name].GetTitle()
        try:
            ja, jb = json.loads(sa), json.loads(sb)
        except ValueError:
            if sa != sb:
                print("TNAMED      :", name, "(title differs, not JSON)"); diff += 1
            continue
        if ja != jb:
            if isinstance(ja, dict) and isinstance(jb, dict):
                for key in sorted(set(ja) | set(jb)):
                    if ja.get(key) != jb.get(key):
                        print("JSON        :", name, key); diff += 1
            else:
                print("JSON        :", name); diff += 1

    print("TOTAL DIFFERENCES:", diff)
    return 1 if diff else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))

#!/usr/bin/env python3
"""Check that the ROOT boxes of a set of certify logs tile one box exactly.

All roots must agree in coordinates 2..dim, and their first-coordinate intervals must chain with no gap
(hi of one == lo of the next, compared as exact floats).  Prints the union box.
With --ref R (the log of an unsliced run, e.g. max_boxes 0), also require the union to equal R's root box.
With --manifest M (written by cloud_instances/aws/split_phi65.py, repeatable): the logs slice_<tag>_s<path>.log of a
split piece <tag> stand for the piece itself.  Its nodes and leaves must form a full binary tree (every node has two
children, leaves at any depth: a hard child may be split again), every leaf must
have exactly one log, and each log's root must equal the box replayed from the piece's root along the leaf's path
(halving the node's coordinate at mid = 0.5*(lo+hi), rounded to nearest, as certify and check_phi bisect).  The
children then tile the piece, which takes their place in the chain (a log of the piece itself as well would overlap
and break the chain).  A child log without its manifest is rejected.
Usage: cover_check.py [--ref R] [--manifest M]... logs...
"""
import json, os, re, sys


def read_root(fn):
    with open(fn) as f:
        head = f.readline().split()
        assert head[0] == "ROOT", fn
        dim = int(head[2].split("=")[1])
        return " ".join(head[1:]), [tuple(float.fromhex(x) for x in f.readline().split()) for _ in range(dim)]


args = sys.argv[1:]
ref, manifests = None, []
while args[:1] in (["--ref"], ["--manifest"]):
    if args[0] == "--ref": ref = read_root(args[1])
    else: manifests.append(json.load(open(args[1])))
    args = args[2:]
logs = {}
for fn in args:
    hdr, box = read_root(fn)
    logs[os.path.basename(fn)] = (hdr, box, fn)

# split pieces: verify each manifest against its child logs, then stand the piece in for them
pieces = []
for m in manifests:
    tag = m["parent"].replace("/", "_")
    root = [(float.fromhex(lo), float.fromhex(hi)) for lo, hi in m["root"]]
    nodes, leaves = m["nodes"], set(m["leaves"])
    # full binary tree: every node's two children are nodes or leaves, and every leaf hangs off a node
    for p in nodes:
        for c in (p + "0", p + "1"):
            if (c in nodes) == (c in leaves): sys.exit("manifest %s: %r is not exactly one of node/leaf" % (tag, c))
    for p in leaves:
        if p and p[:-1] not in nodes: sys.exit("manifest %s: leaf %r has no parent node" % (tag, p))
    if not nodes and leaves != {""}: sys.exit("manifest %s: no nodes" % tag)
    hdrs = set()
    for p in sorted(leaves):
        name = "slice_%s_s%s.log" % (tag, p)
        if name not in logs: sys.exit("split piece %s: no log for leaf %s (gap)" % (tag, p))
        box = list(root)
        for k in range(len(p)):
            t = nodes[p[:k]]
            lo, hi = box[t]
            mid = 0.5 * (lo + hi)
            box[t] = (lo, mid) if p[k] == "0" else (mid, hi)
        if logs[name][1] != box: sys.exit("split piece %s: leaf %s root differs from the replayed box" % (tag, p))
        hdrs.add(logs[name][0])
        del logs[name]
    if len(hdrs) != 1: sys.exit("split piece %s: children have different headers: %s" % (tag, sorted(hdrs)))
    pieces.append((hdrs.pop(), root, "%s (%d sub-logs)" % (m["parent"], len(leaves))))
stray = [n for n in logs if re.search(r"_s[01]*\.log$", n)]
if stray:
    sys.exit("sub-logs without a manifest: %s" % stray[:5])

roots = [(box[0][0], box[0][1], tuple(box[1:]), hdr, fn) for hdr, box, fn in logs.values()]
roots += [(box[0][0], box[0][1], tuple(box[1:]), hdr, name) for hdr, box, name in pieces]
rest, hdr = roots[0][2], roots[0][3]
bad = [r[4] for r in roots if r[2] != rest or r[3] != hdr]
if bad:
    sys.exit("roots disagree outside coordinate 1 (or in n/phi/c): %s" % bad[:5])
roots.sort()
gaps = [(a[4], b[4]) for a, b in zip(roots, roots[1:]) if a[1] != b[0]]
if gaps:
    sys.exit("first-coordinate intervals do not chain: %s" % gaps[:5])
print("%d roots tile y1 in [%s, %s]  (%s)%s" % (len(roots), roots[0][0].hex(), roots[-1][1].hex(), hdr,
      "; %d of them split pieces, each tiled by its sub-logs" % len(pieces) if pieces else ""))
for t, (lo, hi) in enumerate(rest, 2):
    print("  coordinate %d: [%s, %s]" % (t, lo.hex(), hi.hex()))
if ref:
    union = [(roots[0][0], roots[-1][1])] + list(rest)
    if ref[0] != hdr or ref[1] != union:
        sys.exit("union differs from the reference root box (%s)" % ref[0])
    print("union equals the reference root box")

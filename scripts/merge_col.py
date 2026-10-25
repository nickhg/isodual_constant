"""Merge sibling stage-1 collected boxes bottom-up: two boxes equal in 8 coordinates and abutting with equal width in
the ninth are replaced by their union (exactly a box, so coverage is unchanged), repeated to a fixed point.
Usage: python3 merge_col.py [--trees TREES] COLFILE... > merged.col
COLFILEs may be xz-compressed (.xz).  Output lines: the merged box (hex lo hi x 9, the COLLECT format) followed by
"# k" = number of original boxes.  With --trees, TREES gets one line per merged box (same order): its merge tree, a leaf
being "i:l" (COLFILE number i, 0-based, line l) and a merge "(t A B)" (A and B abut in coordinate t, A below), so that
merge_check.py can re-verify every merge from the leaves.  Prints to stderr the count of merged boxes by k."""
import collections, lzma, sys

args = sys.argv[1:]
trees = None
if args and args[0] == "--trees": trees, args = args[1], args[2:]
boxes = []        # (box, k, tree)
for fi, fn in enumerate(args):
    f = lzma.open(fn, "rt") if fn.endswith(".xz") else open(fn)
    for ln, l in enumerate(f):
        v = [float.fromhex(t) for t in l.split()]
        boxes.append((tuple((v[2 * t], v[2 * t + 1]) for t in range(9)), 1, "%d:%d" % (fi, ln)))
changed = True
while changed:
    changed = False
    for t in range(9):
        key = collections.defaultdict(list)
        for i, (b, k, _) in enumerate(boxes):
            key[(b[:t] + b[t + 1:], b[t][1] - b[t][0])].append(i)
        dead, new = set(), []
        for idx in key.values():
            if len(idx) < 2: continue
            bylo = {boxes[i][0][t][0]: i for i in idx}
            for i in sorted(idx, key=lambda i: boxes[i][0][t][0]):
                if i in dead: continue
                j = bylo.get(boxes[i][0][t][1])
                if j is None or j in dead or j == i: continue
                bi, bj = boxes[i][0], boxes[j][0]
                new.append((bi[:t] + ((bi[t][0], bj[t][1]),) + bi[t + 1:], boxes[i][1] + boxes[j][1],
                            "(%d %s %s)" % (t, boxes[i][2], boxes[j][2])))
                dead |= {i, j}
        if new:
            boxes = [x for i, x in enumerate(boxes) if i not in dead] + new
            changed = True
tf = open(trees, "w") if trees else None
for b, k, tr in boxes:
    print(" ".join("%s %s" % (lo.hex(), hi.hex()) for lo, hi in b), "# %d" % k)
    if tf: tf.write(tr + "\n")
c = collections.Counter(k for _, k, _ in boxes)
sys.stderr.write("%d boxes -> %d merged; by size: %s\n" % (sum(k for _, k, _ in boxes), len(boxes), sorted(c.items())))

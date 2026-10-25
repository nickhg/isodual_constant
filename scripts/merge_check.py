"""Independent check of merge_col.py's output: every merged box is the exact union of stage-1 C leaves, and every leaf
is used exactly once.  It re-reads the leaves from the COLFILEs (the STAGE2 output of check_phi on the verified logs),
rebuilds each merged box from its tree, and checks every merge: the two parts are equal in the 8 other coordinates
and abut exactly (hi of A == lo of B) in coordinate t, so their union is a box.  Exact: hex floats, no arithmetic.
Usage: python3 merge_check.py MERGED TREES COLFILE...      (the same COLFILEs, in the same order, as merge_col.py)"""
import lzma, re, sys

merged, trees, cols = sys.argv[1], sys.argv[2], sys.argv[3:]
leaves = []
for fn in cols:
    f = lzma.open(fn, "rt") if fn.endswith(".xz") else open(fn)
    leaves.append([tuple(tuple(float.fromhex(x) for x in l.split()[2 * t:2 * t + 2]) for t in range(9)) for l in f])
used = set()
tok = re.compile(r"\(|\)|[0-9]+:[0-9]+|[0-9]+")


def build(it):
    t = next(it)
    if t == "(":
        c = int(next(it)); a = build(it); b = build(it)
        if next(it) != ")": raise ValueError("bad tree")
        if a[:c] + a[c + 1:] != b[:c] + b[c + 1:]: raise ValueError("parts differ outside coordinate %d" % c)
        if a[c][1] != b[c][0]: raise ValueError("parts do not abut in coordinate %d" % c)
        return a[:c] + ((a[c][0], b[c][1]),) + a[c + 1:]
    i, l = map(int, t.split(":"))
    if (i, l) in used: raise ValueError("leaf %s used twice" % t)
    used.add((i, l))
    return leaves[i][l]


n = 0
for ml, tl in zip(open(merged), open(trees)):
    v = ml.split("#")[0].split()
    box = tuple((float.fromhex(v[2 * t]), float.fromhex(v[2 * t + 1])) for t in range(9))
    if build(iter(tok.findall(tl))) != box: raise SystemExit("merged box %d is not its tree's union" % n)
    n += 1
total = sum(len(L) for L in leaves)
if len(used) != total: raise SystemExit("%d of %d leaves not covered by a merged box" % (total - len(used), total))
print("merge check OK: %d merged boxes = exact unions of all %d C leaves, each used once" % (n, total))

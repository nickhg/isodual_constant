"""Coverage of stage 2 (eta_5 = 6) for a set of stage-1 pieces: every merged box of every piece has a stage-2 log that
check_lift accepted, for exactly that box, or (run_s2.sh's fallback) its two bisected halves are covered in the same way.  For each piece, the chunk files (in job order) must concatenate to
merged/<piece>.col (prep.sh); for each job, its verified .out must report every line i of the chunk as
"logs/<job>/i.log: ver=V check_lift ver=V X <box> OK" with <box> equal (exact hex floats) to that line, and end with
"FAILURES 0".  Together with merge_check.py (the merged boxes are exact unions of all C leaves) and stage 1's coverage,
this covers the pieces.
Usage: python3 cover_s2.py S2DIR VERIFIED VER PIECE...     (S2DIR holds jobs.txt, merged/, chunks/)
VER is the producer's version, or a comma-separated list of accepted producer versions (a producer fixed mid-run);
the checker's may be one of those or one listed in env CHECKVERS (comma-separated), for jobs re-checked with a later
check_lift.  A VERIFIED/<job>.out.recheck, if present, replaces <job>.out.
A VERIFIED/<job>.fix.out (the fix-up, cloud_instances/aws/fix_s2.py: boxes that run_s2.sh's depth-6 bisection left
unproved, and jobs proved box by box) adds its accepted logs to the job's; it must end "FAILURES 0" too, and a job
may have only a .fix.out.  Bisection is accepted to depth MAXDEPTH = 30 (fix_s2.py's limit; run_s2.sh stops at 6)."""
import os, re, sys

s2, ver_dir, ver, pieces = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
pvers = set(ver.split(","))
cvers = pvers | set(filter(None, os.environ.get("CHECKVERS", "").split(",")))
jobs = [l.split() for l in open(os.path.join(s2, "jobs.txt")) if l.strip()]


MAXDEPTH = 30


def box(fields): return tuple(float.fromhex(x) for x in fields)


def halves(b):
    """run_s2.sh's bisection: widest coordinate (ties to the lowest index) at the float midpoint; the two halves are
    equal to b outside coordinate t and meet at m, so their union is exactly b"""
    lo, hi = b[0::2], b[1::2]
    t = max(range(9), key=lambda i: (hi[i] - lo[i], -i)); m = (lo[t] + hi[t]) / 2
    if not lo[t] < m < hi[t]: return None
    out = []
    for a, c in ((lo[t], m), (m, hi[t])):
        l, h = list(lo), list(hi); l[t], h[t] = a, c
        out.append(tuple(x for pair in zip(l, h) for x in pair))
    return out


def covered(rep, ident, b, depth=0):
    """the box b is proved: by an accepted log `ident` for exactly b, or (depth < MAXDEPTH) by its two halves"""
    if ident in rep: return "ok" if rep[ident][2] == b else None
    hv = halves(b) if depth < MAXDEPTH else None
    if hv is None: return None
    return "ok" if all(covered(rep, "%s.%d" % (ident, h), hv[h], depth + 1) == "ok" for h in (0, 1)) else None


bad = 0
for p in pieces:
    pj = [(j, c) for j, c in jobs if re.fullmatch(re.escape(p) + r"_c[0-9]+", j)]
    merged = [l for l in open(os.path.join(s2, "merged", p + ".col")) if l.strip()]
    cat, nb = [], 0
    for j, c in pj:
        lines = [l for l in open(os.path.join(s2, "chunks", c)) if l.strip()]
        cat += lines
        o = os.path.join(ver_dir, j + ".out")
        outs = [o + ".recheck" if os.path.exists(o + ".recheck") else o, os.path.join(ver_dir, j + ".fix.out")]
        outs = [f for f in outs if os.path.exists(f)] or [o]
        rep = {}
        for f in outs:
            out = open(f).read().split("\n")
            for l in out:
                m = re.fullmatch(r"logs/%s/([0-9.]+)\.log: ver=(\S+) check_lift ver=(\S+) X ((?:\S+ ){18})OK" % re.escape(j), l)
                if m:
                    if m.group(1) in rep: print("%s log %s: in two outputs" % (j, m.group(1))); bad += 1
                    rep[m.group(1)] = (m.group(2), m.group(3), box(m.group(4).split()))
            if not any(re.search(r"FAILURES 0$", l) for l in out):
                print("%s: no FAILURES 0 in %s" % (j, f)); bad += 1
        for k, (a, b, _) in rep.items():
            if a not in pvers or b not in cvers: print("%s log %s: version %s/%s" % (j, k, a, b)); bad += 1
        for i, l in enumerate(lines):
            if covered(rep, str(i), box(l.split("#")[0].split())) != "ok":
                print("%s line %d: not covered by accepted logs (whole box or bisected halves)" % (j, i)); bad += 1
        if any(int(k.split(".")[0]) >= len(lines) for k in rep): print("%s: logs beyond the chunk" % j); bad += 1
        nb += len(lines)
    if cat != merged: print("%s: chunks do not concatenate to merged/%s.col" % (p, p)); bad += 1
    print("%s: %d jobs, %d merged boxes" % (p, len(pj), nb))
if bad: raise SystemExit("stage-2 coverage FAILED: %d problems" % bad)
print("stage-2 coverage OK: %d pieces" % len(pieces))

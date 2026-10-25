"""Instance side of the stage-2 fix-up (eta_5 = 6): prove the boxes of FIXES (lines "job id lo_1 hi_1 .. lo_9 hi_9",
from runs/fix_s2_list.py) in parallel, one box per core.  Each box is run with stage2 (TCAP s, default 1800; out of time:
once more with BEAM=1, see prove()); a box
whose log check_lift rejects is re-run with MARGIN 1e-4, then 1e-3; a box out of time, or still rejected, is bisected
(cover_s2.py's rule: widest coordinate, ties to the lowest index, at the float midpoint) and its halves <id>.0, <id>.1
queued, to depth MAXDEPTH (30; run_s2.sh stopped at 6, which left one box unproved).  Box ids are global depths: a gap
box "28.0.0.1" of the stage-2 tree continues as "28.0.0.1.0", ...
Per job, when all its boxes are proved: one check_lift run over all its new logs (sort -V, as run_s2.sh) gives
done/<job>.fix.out; the logs are packed as done/<job>.fix.logs.tar.xz and uploaded to S3 (file S3: s3://bucket/prefix)
with S3's SHA-256 checksum, then done/<job>.fix.sha256 is written (archive, .txt, .out).  A box that fails at
MAXDEPTH stops its job (failed/<job>.out).  Resumable per job.  Progress: progress.log.
Builds: stage2 as run_s2.sh (GITREV = file VERSION); check_lift static as replay_s2.sh (GITREV = file CHECKREV), so
that it is the same binary as the final replay's (sha256 in done/BINARY).
Usage: python3 fix_s2.py FIXES [P]"""
import os, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor

FIXES = sys.argv[1]; P = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
TCAP = os.environ.get("TCAP", "1800"); MAXDEPTH = 30
for d in ("work", "logs", "done", "failed"): os.makedirs(d, exist_ok=True)
REV, CREV, S3 = (open(f).read().strip() for f in ("VERSION", "CHECKREV", "S3"))
sh = lambda c: subprocess.run(c, shell=True, check=True)
sh('gcc -O2 -DGITREV="\\"%s\\"" -o stage2 stage2.c -lglpk -lm' % REV)
sh('gcc -O2 -static -frounding-math -DGITREV="\\"%s\\"" -o check_lift check_lift.c -lm' % CREV)
sh("sha256sum check_lift > done/BINARY; gcc --version | head -1 >> done/BINARY")
lock = threading.Lock()


def log(msg):
    with lock:
        with open("progress.log", "a") as f: f.write("%s %s\n" % (time.strftime("%H:%M:%S", time.gmtime()), msg))


def halves(b):
    lo, hi = b[0::2], b[1::2]
    t = max(range(9), key=lambda i: (hi[i] - lo[i], -i)); m = (lo[t] + hi[t]) / 2
    if not lo[t] < m < hi[t]: return None
    out = []
    for a, c in ((lo[t], m), (m, hi[t])):
        l, h = list(lo), list(hi); l[t], h[t] = a, c
        out.append(tuple(x for pair in zip(l, h) for x in pair))
    return out


def prove(job, ident, b):
    """True if logs/<job>/<ident>.log is a checked certificate of b; else (time out or rejected) False"""
    w = "work/%s/%s" % (job, ident); os.makedirs(w, exist_ok=True)
    with open(w + "/in.col", "w") as f: f.write(" ".join(x.hex() for x in b) + "\n")
    dst = "logs/%s/%s.log" % (job, ident)
    margins = [None, "1e-4", "1e-3"]; k = 0; beam = False
    while k < len(margins):
        env = dict(os.environ, LOGDIR=w, TCAP=TCAP, BOXCAP="100000000")
        if margins[k]: env["MARGIN"] = margins[k]
        if beam: env["BEAM"] = "1"
        r = subprocess.run(["./stage2", w + "/in.col"], env=env, capture_output=True, text=True).stdout
        res = [l for l in r.splitlines() if l.startswith("0 ")][0]
        with lock:
            with open("work/%s.txt" % job, "a") as f: f.write(ident + res[1:] + (" BEAM" if beam else "") + "\n")
        if res.split()[3] != "OK":
            os.remove(w + "/0.log")
            # out of time: once more with stage2's beam search over short-vector flags (BEAM, off by default for
            # speed).  The descent from basis permutations can stall in a local minimum whose eta crosses 6 along a
            # surface, and the producer then refines along that surface without end (the hard paths of 465_480_3_8
            # c004 box 28 and 466_480_{3,6}_8_c000 box 36: exact eta there ~3.9, BEAM proves them in ~90 s)
            if not beam: beam = True; continue
            return False
        os.replace(w + "/0.log", dst)
        c = subprocess.run(["./check_lift", dst], capture_output=True, text=True).stdout
        if c.rstrip().endswith("FAILURES 0"): return True
        os.remove(dst); k += 1
    return False


jobs, todo, state = [], {}, {}
for l in open(FIXES):
    f = l.split()
    if not f: continue
    if f[0] not in todo: jobs.append(f[0]); todo[f[0]] = []
    todo[f[0]].append((f[1], tuple(float.fromhex(x) for x in f[2:20])))
pending = {j: 0 for j in jobs}; failed = set()
pool = ThreadPoolExecutor(P); cv = threading.Condition()


def finish(job):
    out = "done/%s.fix" % job
    sh("./check_lift $(ls logs/%s/*.log | sort -V) > %s.out 2>&1" % (job, out))
    if not open(out + ".out").read().rstrip().endswith("FAILURES 0"):
        sh("cp %s.out failed/%s.out" % (out, job)); log("%s: final check FAILED" % job); return
    sh("cp work/%s.txt %s.txt" % (job, out))
    sh("tar -cJf %s.logs.tar.xz.tmp -C logs/%s . && xz -t %s.logs.tar.xz.tmp && mv %s.logs.tar.xz.tmp %s.logs.tar.xz"
       % (out, job, out, out, out))
    b = S3[5:].split("/", 1)
    for t in range(5):
        if subprocess.run(["aws", "s3api", "put-object", "--bucket", b[0], "--key", "%s/%s.fix.logs.tar.xz" % (b[1], job),
                           "--body", out + ".logs.tar.xz", "--checksum-algorithm", "SHA256"],
                          capture_output=True).returncode == 0: break
        time.sleep(30 * (t + 1))
    else:
        log("%s: S3 upload failed" % job); sh("echo 'S3 upload failed (FAILURES 1)' > failed/%s.out" % job); return
    sh("cd done && sha256sum %s.fix.logs.tar.xz %s.fix.txt %s.fix.out > %s.fix.sha256" % (job, job, job, job))
    log("%s: done (%s)" % (job, open(out + ".out").read().rstrip().splitlines()[-1]))


def task(job, ident, b):
    try:
        ok = prove(job, ident, b)
        if not ok:
            depth = ident.count(".")
            hv = halves(b) if depth < MAXDEPTH else None
            if hv is None:
                log("%s %s: unproved at depth %d" % (job, ident, depth)); failed.add(job)
            else:
                log("%s %s: split (depth %d)" % (job, ident, depth + 1))
                for h in (0, 1): submit(job, "%s.%d" % (ident, h), hv[h])
        else:
            log("%s %s: OK" % (job, ident))
    except Exception as e:
        log("%s %s: error %r" % (job, ident, e)); failed.add(job)
    with cv:
        pending[job] -= 1; last = pending[job] == 0; cv.notify_all()
    if last:
        if job in failed: sh("echo 'unproved box (FAILURES 1)' > failed/%s.out" % job)
        else: finish(job)


def submit(job, ident, b):
    with cv: pending[job] += 1
    pool.submit(task, job, ident, b)


for j in jobs:
    if os.path.exists("done/%s.fix.sha256" % j): continue
    os.makedirs("logs/" + j, exist_ok=True)
    for ident, b in todo[j]: submit(j, ident, b)
with cv:
    cv.wait_for(lambda: all(v == 0 for v in pending.values()))
pool.shutdown(wait=True)
with open("done/FINISHED", "w") as f:
    f.write("finished %s done %d failed %d\n" % (time.strftime("%FT%TZ", time.gmtime()),
            sum(os.path.exists("done/%s.fix.sha256" % j) for j in jobs), len(failed)))

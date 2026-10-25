"""Instance-side driver for check_lift2 on stage-2 archives in S3 (streamed into memory, nothing unpacked to disk).
JOBS lines: "name fraction" with name s2_<job> (eta6_s2/full/<job>.logs.tar.xz) or s2f_<job> (eta6_s2/fix/<job>.fix.logs.tar.xz).
In each archive every (B) leaf is checked, and each other leaf with probability `fraction` (random.Random seeded by the
archive and log names, so the sample is reproducible).  done/<name>.json: counts checked/failed by type and the failures.
Usage: python3 cl2_driver.py JOBS BUCKET [P]"""
import sys, os, io, json, tarfile, random, subprocess, multiprocessing as mp
import check_lift2 as C

def run(job):
    name, frac = job
    out = "done/%s.json" % name
    if os.path.exists(out): return
    j = name.split("_", 1)[1]
    key = ("eta6_s2/fix/%s.fix.logs.tar.xz" % j) if name.startswith("s2f_") else ("eta6_s2/full/%s.logs.tar.xz" % j)
    data = None
    for _ in range(5):
        p = subprocess.run(["aws", "s3", "cp", "s3://%s/%s" % (BUCKET, key), "-"], capture_output=True)
        if p.returncode == 0: data = p.stdout; break
    res = {"name": name, "fraction": frac, "checked": {}, "failed": {}, "failures": [], "error": None}
    if data is None:
        res["error"] = "download failed"
    else:
        try:
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:xz") as t:
                for m in t.getmembers():
                    if not (m.isfile() and m.name.endswith(".log")): continue
                    lines = t.extractfile(m).read().decode().splitlines()
                    rng = random.Random("%s/%s" % (name, m.name))
                    sel = lambda i, typ: typ == "B" or rng.random() < frac
                    for i, (lo, hi, leaf) in enumerate(C.leaves(lines)):
                        if not sel(i, leaf[0]): continue
                        ok = C.prove(C.Cell(lo, hi), leaf); typ = leaf[0]
                        res["checked"][typ] = res["checked"].get(typ, 0) + 1
                        if not ok:
                            res["failed"][typ] = res["failed"].get(typ, 0) + 1; res["failures"].append([m.name, i, typ])
        except Exception as e:
            res["error"] = repr(e)
    json.dump(res, open(out + ".tmp", "w")); os.replace(out + ".tmp", out)

if __name__ == "__main__":
    JOBS, BUCKET = sys.argv[1], sys.argv[2]; P = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    os.makedirs("done", exist_ok=True)
    jobs = [(l.split()[0], float(l.split()[1])) for l in open(JOBS) if l.strip()]
    with mp.Pool(P) as pool: list(pool.imap_unordered(run, jobs))
    open("done/FINISHED", "w").write("finished\n")

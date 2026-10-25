"""The check_lift2 sample (fixed before the first run): every stage-2 job whose archive contains a (B) leaf (the 304 in
runs/replay2/b_jobs.txt), plus 300 of the other jobs drawn with random.Random(20261025); within each archive
cl2_driver.py checks every (B) leaf and each other leaf with probability 0.005, seeded by archive and log name.
Prints the job list ("name fraction"), in the order the first run split it.  Usage: python3 make_sample.py"""
import random
B = [l.split()[0] for l in open("runs/replay2/b_jobs.txt") if l.strip()]
alljobs = ["s2_" + l.split()[0] for l in open("runs/eta6_s2full/jobs.txt") if l.strip()]
alljobs = [j for j in alljobs if j not in set(B)]
rng = random.Random(20261025); R300 = rng.sample(alljobs, 300)
jobs = [(j, 0.005) for j in B] + [(j, 0.005) for j in R300]
rng.shuffle(jobs)
for k in range(2):
    for j in jobs[k::2]: print("%s %g" % j)

P = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
G = Genus(P); reps = G.representatives()
auts = [QuadraticForm(ZZ, 2*matrix(ZZ, R)).number_of_automorphisms() for R in reps]
print("|Aut| of reps:", auts, " sum 1/|Aut| =", sum(1/QQ(a) for a in auts), " mass of genus =", G.mass())
print("minima:", [QuadraticForm(ZZ, 2*matrix(ZZ, R)).theta_series(40).exponents()[1] for R in reps])

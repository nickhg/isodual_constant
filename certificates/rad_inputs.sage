r"""
Inputs of rad_produce.py / rad_check.py (exact radius of the D5 gradient polytope, App. A):
  cf_orth.json   : the 720 matrices C_f = sum_{k<=4} P_{V_k} - 2I in orthonormal coordinates (producer input);
  d5_bases.json  : the 720 Psi-minimal bases of D5 in root coordinates (checker input; the checker rebuilds C_f);
  cap_witness.json: primal witnesses of the 26 cross-polytope LPs of the cap lemma at J - I (checker input).
Also prints the cap-lemma bound s and the cap angle.  Usage: sage rad_inputs.sage
"""
# C_f = sum_{k<=4} P_{V_k} - 2I in orthonormal coordinates of R^5 (D5 in Z^5), for the 720 Psi-minimal flags.
import json
load("psi.sage")
G0 = root_gram('D', 5)
B0 = matrix(QQ, [[1,-1,0,0,0],[0,1,-1,0,0],[0,0,1,-1,0],[0,0,0,1,-1],[0,0,0,1,1]])
assert B0 * B0.T == G0
p0 = psi2_min(G0)[0]; bases = psi2_all(G0, p0)
Cs = set()
for _, U in bases:
    V = U * B0
    C = -2 * identity_matrix(QQ, 5)
    for k in range(1, 5):
        Vk = V.matrix_from_rows(range(k)); C += Vk.T * (Vk * Vk.T).inverse() * Vk
    assert (C * C).trace() == 10 and C.trace() == 0
    Cs.add(tuple(C.list()))
Cs = sorted(Cs); print("distinct C_f:", len(Cs))
J = matrix(QQ, 5, 5, [1] * 25) - identity_matrix(QQ, 5)
vals = [(matrix(QQ, 5, 5, c) * J).trace() for c in Cs]
print("max <C_f, J-I> =", max(vals), "count", vals.count(max(vals)), " |J-I|^2 =", (J * J).trace())
json.dump([[str(x) for x in c] for c in Cs], open("cf_orth.json", "w"))

# ---- cap lemma
Cs = [matrix(QQ, 5, 5, [QQ(x) for x in c]) for c in json.load(open("cf_orth.json"))]
J = matrix(QQ, 5, 5, [1] * 25) - identity_matrix(QQ, 5)
act = [C for C in Cs if (C * J).trace() == 8/3]; assert len(act) == 120
ip = lambda X, Y: (X * Y).trace()
# basis of traceless symmetric orthogonal to J: Gram-Schmidt (Frobenius) of E_ii - E_55, E_ij + E_ji
gens = []
for i in range(4):
    X = matrix(QQ, 5, 5); X[i, i] = 1; X[4, 4] = -1; gens.append(X)
for i in range(5):
    for j in range(i + 1, 5):
        X = matrix(QQ, 5, 5); X[i, j] = X[j, i] = 1; gens.append(X)
B = []
for X in [J] + gens:
    Y = X
    for Z in B: Y = Y - ip(Y, Z) / ip(Z, Z) * Z
    if Y != 0: B.append(Y)
assert len(B) == 14; T = B[1:]                 # tangent basis, orthogonal, |T_k|^2 rational
G = [vector(QQ, [ip(C, Tk) for Tk in T]) for C in act]     # coordinates <C_f, T_k>
def axis_t(d):   # largest t with t*d in conv(G), d in coefficient space of the T_k
    P = MixedIntegerLinearProgram(maximization=True, solver="PPL")
    lam = P.new_variable(nonnegative=True); t = P.new_variable()
    P.add_constraint(sum(lam[f] for f in range(len(G))) == 1)
    for k in range(13): P.add_constraint(sum(lam[f] * G[f][k] for f in range(len(G))) == t[0] * d[k])
    P.set_objective(t[0]); return QQ(P.solve())
# functional g = sum_k y_k T_k^* ... in coefficient coordinates y_k = <C, T_k>, the functional norm is
# |g|^2 = sum_k y_k^2 / |T_k|^2  (T orthogonal).  Axis e_k: unit vector in y-coordinates; its norm 1/|T_k|.
inv2 = 0
for k in range(13):
    e = vector(QQ, [1 if i == k else 0 for i in range(13)])
    tp, tm = axis_t(e), axis_t(-e)
    assert tp > 0 and tm > 0
    a2 = min(tp, tm) ** 2 / ip(T[k], T[k])
    inv2 += 1 / a2
    print("axis %2d  t+ %s  t- %s  a %.4f" % (k, tp, tm, float(sqrt(a2))), flush=True)
s2 = 1 / inv2; s = sqrt(RIF(s2)); r0 = RIF(4) / (3 * sqrt(RIF(5)))
print("s^2 >=", s2, " s >=", s, " r0 =", r0, "  s/r0 >=", s / r0)
th = 2 * arctan(s / r0); print("theta0 =", th * 180 / RIF(pi), "deg")

# ---- witnesses and bases
G0 = root_gram("D", 5); p0 = psi2_min(G0)[0]; bases = psi2_all(G0, p0)
json.dump([[int(x) for x in U.list()] for _, U in bases], open("d5_bases.json", "w"))

wit = []
for k in range(13):
    for sgn in (1, -1):
        d = vector(QQ, [sgn if i == k else 0 for i in range(13)])
        P = MixedIntegerLinearProgram(maximization=True, solver="PPL")
        lam = P.new_variable(nonnegative=True); t = P.new_variable()
        P.add_constraint(sum(lam[f] for f in range(len(G))) == 1)
        for kk in range(13): P.add_constraint(sum(lam[f] * G[f][kk] for f in range(len(G))) == t[0] * d[kk])
        P.set_objective(t[0]); tv = QQ(P.solve()); lv = P.get_values(lam)
        sup = {f: QQ(lv[f]) for f in range(len(G)) if QQ(lv[f]) != 0}
        assert sum(sup.values()) == 1 and sum(l * G[f] for f, l in sup.items()) == tv * d
        # record the active flag by its C_f matrix (index into the active list order of cf_orth.json)
        wit.append({"k": int(k), "sign": int(sgn), "t": str(tv), "lam": {str(f): str(l) for f, l in sup.items()}})
json.dump({"witnesses": wit, "tan_half_theta0": "27/100"}, open("cap_witness.json", "w"))
print("bases", len(bases), "witnesses", len(wit))

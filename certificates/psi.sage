# psi.sage -- exact minimisation of the Psi-measure (LLL potential) of a lattice.
#
#   Psi(B) = prod_i |b_i*|^(n+1-i) = prod_k det(L_k),   L_k = span(b_1..b_k)
#   eta(L) = Psi_min(L)^2 / Delta^(n+1),                Delta = covolume
#   eta_n  = max_L eta(L)          (ePrint 2007/105: eta_2 = gamma_2 = sqrt(4/3), eta_3 = 3/2)
#
# Lattices are given by a positive-definite rational Gram matrix G, and
# everything is computed in squared lengths, so all arithmetic is exact.
#
# Search: Psi^2(M) = min over primitive v of (v G v)^m * Psi^2(pi_v M), where m =
# dim M and pi_v projects orthogonally to v.  A Psi-minimal basis, after size
# reduction (which leaves Psi unchanged), is LLL-reduced with delta = 1: each
# consecutive 2-block must be Lagrange-reduced or a swap would lower Psi.  So we
# only enumerate LLL(delta=1) bases:
#   |b_1*|^2 <= (4/3)^((m-1)/2) det(M)^(1/m)     (bound on the first vector of each projection)
#   |b_{i+1}*|^2 >= (3/4) |b_i*|^2               (Lovasz with delta=1, |mu| <= 1/2)
# and prune with the Lovasz chain lower bound on the rest of the basis.

from sage.all import QQ, ZZ, matrix, vector, gcd, pari, prod, CartanMatrix


def _complete_basis(v):
    """Unimodular integer matrix whose first row is the primitive vector v."""
    m = len(v)
    D, S, T = matrix(ZZ, [list(v)]).smith_form()   # D = S * v * T = (1, 0, ..., 0)
    W = T.inverse()                                 # v * T = S^-1 * e1, so row 0 of W is +-v
    if vector(ZZ, W.row(0)) != vector(ZZ, v):
        W.rescale_row(0, -1)
    assert W.det() in (1, -1) and vector(ZZ, W.row(0)) == vector(ZZ, v)
    return W


def _project(G, v):
    """Gram matrix of the projection of the lattice orthogonal to v (dimension m-1)."""
    W = _complete_basis(v)
    H = W * G * W.transpose()
    a = H[0, 0]
    c = H[1:, 0]
    return H[1:, 1:] - (c * c.transpose()) / a


def _short_primitive(G, bound):
    """Primitive v (one of each +-v) with v G v <= bound, sorted by norm."""
    den = G.denominator()
    Gi = matrix(ZZ, G * den)
    B = (QQ(bound) * den).floor()
    if B <= 0:
        return []
    res = pari(Gi).qfminim(B, None, 2)   # flag 2: exact, no count limit
    out = []
    for col in res[2]:
        v = vector(ZZ, [ZZ(x) for x in col])
        if gcd(list(v)) == 1:
            out.append((v * G * v, v))
    out.sort(key=lambda t: t[0])
    return out


def _rest_lb(first_sq, m):
    """Lower bound on Psi^2 of an m-dim LLL(delta=1) tail whose first |b*|^2 >= first_sq."""
    return prod((QQ(3) / 4) ** (j * (m - j)) * first_sq ** (m - j) for j in range(m))


def _lll_gram(G):
    """Gram matrix of an LLL-reduced basis of G (exact, via pari qflllgram on a scaled copy)."""
    den = G.denominator()
    U = matrix(ZZ, pari(matrix(ZZ, G * den)).qflllgram())
    return U.transpose() * G * U


def _gs_profile(G):
    """Squared Gram-Schmidt norms r_1..r_m of the basis with Gram matrix G."""
    dets = [QQ(1)] + [G[:k, :k].det() for k in range(1, G.nrows() + 1)]
    return [dets[k] / dets[k - 1] for k in range(1, G.nrows() + 1)]


def psi2_lb(G):
    """Lower bound on Psi_min^2 from Howgrave-Graham 1998 Lemma 3.2.5 (= Pataki-Tural Lemma 1).

    Every rank-k sublattice has Gram determinant >= the product of the k smallest r_i of
    any basis, so Psi^2 = prod_k det(L_k) >= prod_k (product of the k smallest r_i).
    We use an LLL-reduced basis to make the bound reasonably tight.
    """
    r = sorted(_gs_profile(_lll_gram(G)))
    out, acc = QQ(1), QQ(1)
    for x in r[:-1]:
        acc *= x
        out *= acc
    return out * G.det()


def psi2_of(G):
    return prod(G[:k, :k].det() for k in range(1, G.nrows() + 1))


def psi2_min(G, prev_sq=0, best=None, stats=None, use_lb=True):
    """(Psi_min^2, flag of squared Gram-Schmidt norms) for Gram matrix G.

    prev_sq: |b*|^2 of the previous basis vector (Lovasz constraint); best: prune bound.
    At top level (best=None) the search is seeded with Psi^2 of an LLL-reduced basis.
    use_lb: prune subtrees with the Lemma 3.2.5 lower bound psi2_lb.
    """
    m = G.nrows()
    if m == 1:
        return G[0, 0], [G[0, 0]]
    if stats is not None:
        stats['nodes'] = stats.get('nodes', 0) + 1
    if best is None and prev_sq == 0 and m >= 3:
        best = psi2_of(_lll_gram(G)) * (1 + QQ(1) / 10 ** 9)
    det = G.det()
    # LLL(delta=1) bound on the first vector, squared: (4/3)^((m-1)/2) det^(1/m);
    # compare via powers to stay exact: x^m <= (4/3)^(m(m-1)/2) det.
    lll_cap = (QQ(4) / 3) ** (m * (m - 1) / 2) * det
    cap = (lll_cap.n(200) ** (1 / QQ(m))) * (1 + 10 ** -30)
    best_val, best_flag = (None, None) if best is None else (best, None)
    lo = QQ(3) / 4 * prev_sq
    for n2, v in _short_primitive(G, cap):
        if n2 ** m > lll_cap or n2 < lo:
            continue
        # prune: n2^m * (tail lower bound with first |b*|^2 >= 3/4 n2)
        if best_val is not None and n2 ** m * _rest_lb(QQ(3) / 4 * n2, m - 1) >= best_val:
            break   # list is sorted by n2 and the bound is increasing in n2
        P = _project(G, v)
        if use_lb and best_val is not None and m >= 3 and n2 ** m * psi2_lb(P) >= best_val:
            if stats is not None:
                stats['lb_pruned'] = stats.get('lb_pruned', 0) + 1
            continue
        sub_best = None if best_val is None else best_val / n2 ** m
        r = psi2_min(P, n2, sub_best, stats, use_lb)
        if r[0] is None:
            continue
        val = n2 ** m * r[0]
        if best_val is None or val < best_val:
            best_val, best_flag = val, [n2] + r[1]
    return best_val, best_flag


def psi2_all(G, cap, prev_sq=0, M=None, out=None):
    """All LLL(delta=1) bases U (rows = integer coords in the input basis) with Psi^2 <= cap.

    Returns a list of (Psi^2, U).  Used to find the active set of near-minimal flags.
    """
    m = G.nrows()
    if M is None:
        M = identity_matrix(ZZ, m)
    if out is None:
        out = []
    if m == 1:
        if G[0, 0] <= cap:
            out.append((G[0, 0], M))
        return out
    det = G.det()
    lll_cap = (QQ(4) / 3) ** (m * (m - 1) / 2) * det
    rad = (lll_cap.n(200) ** (1 / QQ(m))) * (1 + 10 ** -30)
    lo = QQ(3) / 4 * prev_sq
    for n2, v in _short_primitive(G, rad):
        if n2 ** m > lll_cap or n2 < lo:
            continue
        if n2 ** m * _rest_lb(QQ(3) / 4 * n2, m - 1) > cap:
            break
        W = _complete_basis(v)
        P = _project(G, v)
        WM = W * M
        sub = psi2_all(P, cap / n2 ** m, n2, WM.matrix_from_rows(range(1, m)), [])
        for val, U in sub:
            out.append((n2 ** m * val, WM.matrix_from_rows([0]).stack(U)))
    return out


def lambda1(G):
    """Exact squared minimum of the lattice with Gram matrix G."""
    return _short_primitive(G, _lll_gram(G)[0, 0])[0][0]


def psi2_hkz(G):
    """Smallest Psi^2 over HKZ bases (first vector of every projection is a shortest vector)."""
    m = G.nrows()
    if m == 1:
        return G[0, 0]
    sv = _short_primitive(G, lambda1(G))
    lam = sv[0][0]
    return min(lam ** m * psi2_hkz(_project(G, v)) for n2, v in sv if n2 == lam)


def eta(G, stats=None):
    """eta(L) = Psi_min^2 / det(G)^((n+1)/2), with the minimising squared GS profile."""
    n = G.nrows()
    p2, flag = psi2_min(matrix(QQ, G), stats=stats)
    return p2 / matrix(QQ, G).det() ** (QQ(n + 1) / 2), flag


def root_gram(t, n):
    """Gram matrix (min norm 2) of the root lattice of Cartan type t_n, e.g. ('E', 8)."""
    return matrix(QQ, CartanMatrix([t, n]))

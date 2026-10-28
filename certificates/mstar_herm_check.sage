# Check that h(x,y) = (1/3) sum_k <g^k x, y> w^k is Hermitian for omega acting as g, with trace <,>.
P = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
g = matrix(ZZ, [[0,-1,1,0],[1,-1,1,-1],[0,0,0,-1],[0,0,1,-1]])
K.<w> = CyclotomicField(3)
h = lambda x, y: sum((g^k*x)*P*y * w^k for k in range(3)) / 3
E = identity_matrix(ZZ, 4).columns()
ok = all(h(g*x, y) == w^-1*h(x, y) or h(g*x, y) == w*h(x, y) for x in E for y in E)
lin = set()
for x in E:
    for y in E:
        lin.add("w^-1" if h(g*x, y) == w^-1*h(x, y) else ("w" if h(g*x, y) == w*h(x, y) else "?"))
print("scalar on first argument:", lin)
print("h(y,x) = conj h(x,y):", all(h(y, x) == h(x, y).conjugate() for x in E for y in E))
print("trace = <,>:", all(h(x, y).trace() == x*P*y for x in E for y in E))
print("h(x, g y) vs h(x,y):", set(("w" if h(x, g*y) == w*h(x, y) else "w^-1" if h(x, g*y) == w^-1*h(x, y) else "?") for x in E for y in E))

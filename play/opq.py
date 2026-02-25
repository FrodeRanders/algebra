import algebrapy as alg
F = alg.Fp(23)
a = F.elem(9)
print(F.is_quadratic_residue(a), F.sqrt(a))
g = F.elem(5)
h = F.pow(g, 7)
print(F.dlog(g, h))  # 7

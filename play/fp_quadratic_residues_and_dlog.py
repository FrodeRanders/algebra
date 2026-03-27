import algebrapy as alg

F = alg.Fp(23)
a = F.elem(9)
print("is 9 a quadratic residue?", F.is_quadratic_residue(a))
print("a square root of 9 =", F.sqrt(a))

g = F.elem(5)
h = F.pow(g, 7)
print("g =", g)
print("h = g^7 =", h)
print("discrete log of h base g =", F.dlog(g, h))

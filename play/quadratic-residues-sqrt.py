import algebrapy as alg
F = alg.Fp(23)
for x in [0,1,2,3,4,5,9,10]:
    a = F.elem(x)
    print(x, F.is_quadratic_residue(a), F.sqrt(a))

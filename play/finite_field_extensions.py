import algebrapy as alg

f = alg.PolyFp(2, [1,1,0,1])  # x^3 + x + 1
print(f, f.degree(), f.is_irreducible())

K = alg.Fq(2, [1,1,0,1])
x = K.elem([0,1])             # x
print(K.size(), x)
print([K.pow(x, i).coeffs() for i in range(0, 8)])

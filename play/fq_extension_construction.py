import algebrapy as alg

f = alg.PolyFp(2, [1, 1, 0, 1])  # coefficients low -> high, so this is 1 + x + x^3
print("modulus polynomial =", f)
print("degree =", f.degree())
print("irreducible? =", f.is_irreducible())

K = alg.Fq(2, [1, 1, 0, 1])  # same low -> high convention
x = K.elem([0, 1])
print("field size =", K.size())
print("distinguished element x =", x)
print("powers of x =", [K.pow(x, i).coeffs() for i in range(8)])

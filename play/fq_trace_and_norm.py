import algebrapy as alg

K = alg.Fq(2, [1, 1, 0, 1])  # coefficients low -> high: 1 + x + x^3
a = K.elem([1, 0, 1])
print("a =", a)
print("trace(a) =", a.trace())
print("norm(a)  =", a.norm())

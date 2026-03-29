import algebrapy as alg

K = alg.Fq(3, [2, 0, 1])  # coefficients low -> high: 2 + x^2 over F3
a = K.elem([1, 2])  # 1 + 2x

print("a        =", a)
print("trace(a) =", a.trace())
print("norm(a)  =", a.norm())

gens = K.primitive_elements()
print("primitive elements:", gens)
print("orders:", [g.mul_order() for g in gens])

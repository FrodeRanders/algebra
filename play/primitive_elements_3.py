import algebrapy as alg

K = alg.Fq(3, [2,0,1])      # x^2 + 2 over F3 (monic)
a = K.elem([1,2])           # 1 + 2x

print("a        =", a)
print("trace(a) =", a.trace())
print("norm(a)  =", a.norm())

gens = K.primitive_elements()
print("primitive elements (some):", gens[:10])
#
# Expected: ValueError: failed to find multiplicative order (unexpected)
# since this is a quotient ring and not a field, so the assumption
# "multiplicative group of nonzero elements has order q-1” is false.
#

import algebrapy as alg

K = alg.Fq(2, [1, 1, 0, 1])  # GF(8)
gens = K.primitive_elements()
print("primitive elements:", gens)
print("count:", len(gens))
print("orders:", [g.mul_order() for g in gens])

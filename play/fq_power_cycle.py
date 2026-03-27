import algebrapy as alg

f = alg.PolyFp(2, [1, 1, 0, 1])
print("irreducible?", f.is_irreducible())

K = alg.Fq(2, [1, 1, 0, 1])
x = K.elem([0, 1])
print("field size =", K.size())
for i in range(8):
    print(f"x^{i} =", K.pow(x, i))

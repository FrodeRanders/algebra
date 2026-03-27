import algebrapy as alg

F = alg.Fp(7)
a = F.elem(3)
b = F.elem(5)

print("field modulus =", F.modulus())
print("elements =", F.elements())
print("nonzero elements =", F.nonzero_elements())
print("a =", a, "b =", b)
print("a+b =", F.add(a, b))
print("a*b =", F.mul(a, b))
print("inv(a) =", F.inv(a))
print("a^6 =", F.pow(a, 6))
print("order(a) =", F.mul_order(a))

# Quick check: Fermat's little theorem for all nonzero elements.
for x in F.nonzero_elements():
    assert int(F.pow(x, 6)) == 1
print("Fermat check OK")

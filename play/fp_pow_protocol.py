import algebrapy as alg

F = alg.Fp(7)
a = F.elem(3)

print("a =", a)
print("a ** 6 =", a**6)
print("a ** -1 =", a**-1)
print("pow(a, 6) =", pow(a, 6))

import algebrapy as alg

F = alg.Fp(7)
a = F.elem(3)
b = F.elem(5)

print("a + b =", a + b)
print("a - b =", a - b)
print("a * b =", a * b)
print("-a =", -a)
print("a ** 6 =", a**6)
print("a ** -1 =", a**-1)
print("a / b =", a / b)

print("3 + a =", 3 + a)
print("3 * a =", 3 * a)
print("3 - a =", 3 - a)
print("3 / a =", 3 / a)

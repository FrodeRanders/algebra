import algebrapy as alg

R = alg.Zn(12)
a = R.elem(5)
b = R.elem(6)

print("modulus =", R.modulus())
print("elements =", R.elements())
print("a =", a, "b =", b)
print("a+b =", a + b)
print("a*b =", a * b)
print("a^3 =", a ** 3)
print("units =", R.units())
print("a is unit?", a.is_unit())
print("b is unit?", b.is_unit())
print("inv(a) =", a.inv())

try:
    print("inv(b) =", b.inv())
except Exception as exc:
    print("inv(b) failed:", exc)

print("3 + a =", 3 + a)
print("3 * a =", 3 * a)
print("3 / a =", 3 / a)

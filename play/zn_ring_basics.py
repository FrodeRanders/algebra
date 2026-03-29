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
print("zero divisors =", R.zero_divisors())
print("R is integral domain?", R.is_integral_domain())
print("a is unit?", a.is_unit())
print("b is unit?", b.is_unit())
print("a is zero divisor?", a.is_zero_divisor())
print("b is zero divisor?", b.is_zero_divisor())
print("inv(a) =", a.inv())

try:
    print("inv(b) =", b.inv())
except Exception as exc:
    print("inv(b) failed:", exc)

print("3 + a =", 3 + a)
print("3 * a =", 3 * a)
print("3 / a =", 3 / a)

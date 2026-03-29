import algebrapy as alg

R = alg.Zn(12)
u = R.elem(5)
z = R.elem(6)
d = R.elem(4)
e = R.elem(3)

print("modulus =", R.modulus())
print("elements =", R.elements())
print("u =", u, "z =", z)
print("u+z =", u + z)
print("u*z =", u * z)
print("u^3 =", u ** 3)
print("units =", R.units())
print("zero divisors =", R.zero_divisors())
print("R is integral domain?", R.is_integral_domain())
print("u is unit?", u.is_unit())
print("z is unit?", z.is_unit())
print("u is zero divisor?", u.is_zero_divisor())
print("z is zero divisor?", z.is_zero_divisor())
print("d =", d, "e =", e)
print("d*e =", d * e, "so both are zero divisors")
print("inv(u) =", u.inv())

try:
    print("inv(z) =", z.inv())
except Exception as exc:
    print("inv(z) failed:", exc)

print("3 + u =", 3 + u)
print("3 * u =", 3 * u)
print("3 / u =", 3 / u)

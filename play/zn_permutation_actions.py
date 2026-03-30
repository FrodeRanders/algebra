import algebrapy as alg

R = alg.Zn(12)
u = R.elem(5)
z = R.elem(6)
b = R.elem(3)

print("R =", R)
print("units =", R.units())
print("zero divisors =", R.zero_divisors())
print()

tu = R.add_perm(b)
mu = R.mul_perm(u)
au = R.affine_perm(u, b)

print("translation x -> x + 3")
print("perm =", tu)
print("images =", tu.as_images())
print("cycles =", tu.cycles())
print("cycle type =", tu.cycle_type())
print("order =", tu.order())
print()

print("multiplication x -> 5x")
print("perm =", mu)
print("images =", mu.as_images())
print("cycles =", mu.cycles())
print("cycle type =", mu.cycle_type())
print("order =", mu.order())
print()

print("affine action x -> 5x + 3")
print("perm =", au)
print("images =", au.as_images())
print("cycles =", au.cycles())
print("cycle type =", au.cycle_type())
print("order =", au.order())
print()

print("The unit group acting by multiplication:")
for perm in R.unit_group_perms():
    print(
        perm.as_images(),
        "cycles =",
        perm.cycles(),
        "cycle type =",
        perm.cycle_type(),
        "order =",
        perm.order(),
    )
print()

try:
    print("x -> 6x =", R.mul_perm(z))
except Exception as exc:
    print("x -> 6x failed:", exc)

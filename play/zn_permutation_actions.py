import algebrapy as alg


def show_perm(label, perm):
    print(label)
    print("  perm =", perm)
    print("  images =", perm.as_images())
    print("  cycles =", perm.cycles())
    print("  cycle type =", perm.cycle_type())
    print("  order =", perm.order())
    print()


def show_action(title, ring, generators):
    print(title)
    print("modulus =", ring.modulus())
    print("units =", ring.units())
    print("zero divisors =", ring.zero_divisors())
    sn = alg.Sn(ring.modulus())
    print("orbits =", sn.orbits(generators))
    print()


# Composite modulus: the unit group action is not transitive on nonzero residues.
R12 = alg.Zn(12)
u5 = R12.elem(5)
u7 = R12.elem(7)
z6 = R12.elem(6)
b3 = R12.elem(3)

t12 = R12.add_perm(b3)
m12 = R12.mul_perm(u5)
a12 = R12.affine_perm(u5, b3)

show_perm("Z/12Z translation x -> x + 3", t12)
show_perm("Z/12Z multiplication x -> 5x", m12)
show_perm("Z/12Z affine action x -> 5x + 3", a12)

show_action("Unit-group action on Z/12Z", R12, [R12.mul_perm(u5), R12.mul_perm(u7)])

try:
    print("x -> 6x =", R12.mul_perm(z6))
except Exception as exc:
    print("x -> 6x failed:", exc)
print()

# Prime modulus: multiplication by a primitive element acts transitively on nonzero residues.
R7 = alg.Zn(7)
g = R7.elem(3)  # primitive modulo 7
t7 = R7.add_perm(R7.elem(1))
m7 = R7.mul_perm(g)
a7 = R7.affine_perm(g, R7.elem(1))

show_perm("Z/7Z translation x -> x + 1", t7)
show_perm("Z/7Z multiplication x -> 3x", m7)
show_perm("Z/7Z affine action x -> 3x + 1", a7)

show_action("Cyclic unit-group action on Z/7Z", R7, [m7])

print("orbit of 1 under multiplication by 3 mod 7 =", alg.Sn(7).orbit(1, [m7]))

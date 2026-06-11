import algebrapy as alg

print("=== Polynomial arithmetic over F5 ===")
f = alg.PolyFp(5, [1, 2, 1])   # 1 + 2x + x^2
g = alg.PolyFp(5, [0, 1])       # x
h = alg.PolyFp(5, [2, 0, 1])    # 2 + x^2

print("f =", f)
print("g =", g)
print("h =", h)

print()
print("f + g =", f + g)
print("f - g =", f - g)
print("f * g =", f * g)
print("g * h =", g * h)
print("-f    =", -f)
print("g**3  =", g**3)

print()
print("=== Polynomial division ===")
q, r = h.div_rem(f)
print("h / f : quotient =", q, ", remainder =", r)
print("check: f*q + r =", f * q + r)

print()
print("=== Polynomial gcd and egcd ===")
a = alg.PolyFp(2, [1, 0, 1])   # 1 + x^2  (note: (1+x)^2 = 1+x^2 in F2)
b = alg.PolyFp(2, [1, 1, 0, 1]) # 1 + x + x^3
g, s, t = a.egcd(b)
lhs = s * a + t * b
print("a =", a)
print("b =", b)
print("gcd(a,b) =", g.monic())
print("Bezout: s*a + t*b =", lhs.monic())

print()
print("=== Irreducibility over F2 ===")
for coeffs, name in [
    ([1, 1, 1], "1 + x + x^2"),
    ([1, 0, 1], "1 + x^2"),
    ([1, 1, 0, 1], "1 + x + x^3"),
    ([1, 0, 1, 1], "1 + x^2 + x^3"),
]:
    p = alg.PolyFp(2, coeffs)
    print(f"  {name}: irreducible = {p.is_irreducible()}")

print()
print("=== Polynomial evaluation over F7 ===")
f = alg.PolyFp(7, [1, 2, 3])  # 1 + 2x + 3x^2
print("f(x) = 1 + 2x + 3x^2")
for x in range(7):
    print(f"  f({x}) = {f.eval(x)}")

print()
print("=== Evaluation of a known polynomial ===")
# Lagrange interpretation: f(0)=1, f(1)=3, f(2)=2 in F5
# f(x) = (1 - x)(2 - x) * 3^{-1} + ...  (just checking eval matches)
f = alg.PolyFp(5, [1, 4, 1])  # 1 + 4x + x^2
for x in range(5):
    print(f"  f({x}) = {f.eval(x)}")

import algebrapy as alg

F = alg.Fp(7)
a = F.elem(3)
b = F.elem(5)

print(a + b)        # 1 (mod 7)
print(a * b)        # 1 (mod 7)
print(-a)           # 4 (mod 7)
print(a ** 6)       # 1 (mod 7)
print(a ** -1)      # 5 (mod 7)
print(a / b)        # 2 (mod 7) because 3*5^{-1}=3*3=9=2

print(3 + a)        # 6 (mod 7)
print(3 * a)        # 2 (mod 7)
print(3 - a)        # 0 (mod 7)
print(3 / a)        # 1 (mod 7) because 3*3^{-1}=3*5=15=1

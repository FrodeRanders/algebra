import algebrapy as alg
F = alg.Fp(7)
a = F.elem(3)
print(a**6)      # 1 (mod 7)
print(a**-1)     # 5 (mod 7)
print(pow(a, 6)) # 1 (mod 7)

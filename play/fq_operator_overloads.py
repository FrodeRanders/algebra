import algebrapy as alg

K = alg.Fq(2, [1, 1, 0, 1])  # 1 + x + x^3
a = K.elem([1, 0, 1])  # 1 + x^2
b = K.elem([0, 1])  # x

print("a:", a)
print("b:", b)
print("a+b:", a + b)
print("a-b:", a - b)
print("a*b:", a * b)
print("b**7:", b**7)
print("a/b:", a / b)
print("a**-1:", a**-1)
print("3 + a:", 3 + a)
print("2 * b:", 2 * b)

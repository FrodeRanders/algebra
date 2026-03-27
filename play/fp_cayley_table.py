import algebrapy as alg

F = alg.Fp(5)
elems = [int(e) for e in F.elements()]
mul_tbl = [[int(F.mul(F.elem(i), F.elem(j))) for j in elems] for i in elems]

header = "   | " + " ".join(f"{j:>2}" for j in elems)
rule = "-" * len(header)

print("multiplication table for F_5:")
print(header)
print(rule)
for i, row in zip(elems, mul_tbl):
    print(f"{i:>2} | " + " ".join(f"{value:>2}" for value in row))

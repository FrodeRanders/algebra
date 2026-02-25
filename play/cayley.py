import pandas as pd
import algebrapy as alg

F = alg.Fp(5)
elems = [int(e) for e in F.elements()]
mul_tbl = [[int(F.mul(F.elem(i), F.elem(j))) for j in elems] for i in elems]
df = pd.DataFrame(mul_tbl, index=elems, columns=elems)

print(df)
input("\nPress Enter to exit...")

import algebrapy as alg

S4 = alg.Sn(4)
t = alg.Perm.cycle(4, [0, 1])
c = alg.Perm.cycle(4, [0, 1, 2, 3])
G = S4.generated([t, c])

print("G =", G)
print("order(G) =", G.order())
print("is abelian? =", G.is_abelian())
print()

for p in [2, 3]:
    sylows = G.sylow_p_subgroups(p)
    print(f"Sylow {p}-subgroups:", len(sylows))
    print("orders =", [H.order() for H in sylows])
    print("normal? =", [H.is_normal_in(G) for H in sylows])
    print()

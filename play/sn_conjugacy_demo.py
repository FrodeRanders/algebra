import algebrapy as alg

S4 = alg.Sn(4)
t = alg.Perm.cycle(4, [0, 1])
c = alg.Perm.cycle(4, [0, 1, 2])
g = alg.Perm.cycle(4, [0, 3, 2, 1])

print("S4")
print("t =", t, "cycle notation =", t.cycle_notation())
print("c =", c, "cycle notation =", c.cycle_notation())
print()

print("conjugating t by g")
tg = t.conjugate_by(g)
print("g =", g.cycle_notation())
print("g t g^-1 =", tg.cycle_notation())
print("cycle type preserved?", tg.cycle_type() == t.cycle_type())
print()

print("conjugacy class size of a transposition =", S4.conjugacy_class_size(t, 24))
print("conjugacy class size of a 3-cycle =", S4.conjugacy_class_size(c, 24))

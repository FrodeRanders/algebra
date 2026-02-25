import algebrapy as alg

F = alg.Fp(11)
G = [F.elem(i) for i in range(1, 11)]
orders = {int(a): a.mul_order() for a in G}
print(orders)

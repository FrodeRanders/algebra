import algebrapy as alg

F = alg.Fp(11)
units = [F.elem(i) for i in range(1, 11)]
orders = {int(a): a.mul_order() for a in units}

print("multiplicative orders in F_11^*:")
for value, order in orders.items():
    print(f" {value}: order {order}")

primitive = [value for value, order in orders.items() if order == 10]
print("primitive elements =", primitive)

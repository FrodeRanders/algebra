import algebrapy as alg
f = alg.PolyFp(2, [1,1,0,1])
print("irreducible?", f.is_irreducible())
K = alg.Fq(2, [1,1,0,1])
x = K.elem([0,1])
print("K size", K.size())
print("powers of x:", [K.pow(x,i).coeffs() for i in range(8)])

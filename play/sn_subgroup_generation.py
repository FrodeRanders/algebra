import algebrapy as alg

S3 = alg.Sn(3)
sigma = alg.Perm.cycle(3, [0, 1, 2])
tau = alg.Perm.cycle(3, [0, 1])
H = S3.generated([sigma, tau])
print("generators:", sigma, tau)
print("generated size =", len(H))
print("subgroup =", H)

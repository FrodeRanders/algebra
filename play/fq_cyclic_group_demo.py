import algebrapy as alg


def main():
    f = alg.PolyFp(2, [1, 1, 0, 1])  # coefficients low -> high: 1 + x + x^3
    print("modulus:", f)
    print("irreducible?", f.is_irreducible())

    K = alg.Fq(2, [1, 1, 0, 1])  # GF(8), same low -> high convention
    q = K.size()
    print("field size q =", q)
    print("|Fq*| =", q - 1)

    gens = K.primitive_elements()
    print("primitive elements:", gens)
    if not gens:
        print("no primitive element found")
        return

    g = gens[0]
    print("chosen generator g =", g)
    print("order(g) =", K.mul_order(g))

    powers = [K.pow(g, i) for i in range(q - 1)]
    print("powers g^i for i = 0..q-2:")
    for i, value in enumerate(powers):
        print(f"  g^{i} =", value)

    nonzero = [a for a in K.elements(None) if a]
    print("all nonzero elements:", nonzero)
    print("powers cover all nonzero elements?", set(map(repr, powers)) == set(map(repr, nonzero)))


if __name__ == "__main__":
    main()

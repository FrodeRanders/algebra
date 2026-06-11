import pytest
import algebrapy as alg


def gf8():
    return alg.Fq(2, [1, 1, 0, 1])


def gf9():
    return alg.Fq(3, [1, 1, 0, 0, 1])


def test_fq_construction_requires_monic_irreducible():
    # 1 + x + x^3 is monic and irreducible over F2, so it defines GF(2^3).
    k = alg.Fq(2, [1, 1, 0, 1])
    assert k.p() == 2
    assert k.degree() == 3
    assert k.size() == 8

    # F4 is not a prime base field, so F4[x]/(f) is not this constructor's model.
    with pytest.raises(ValueError):
        alg.Fq(4, [1, 1, 1])
    # 1 + x^2 = (x + 1)^2 over F2, so the quotient would have zero divisors.
    with pytest.raises(ValueError):
        alg.Fq(2, [1, 0, 1])
    # 1 + 2*x^2 is irreducible over F3 but not monic, and reduction assumes monic f.
    with pytest.raises(ValueError):
        alg.Fq(3, [1, 0, 2])


def test_fq_size_and_modulus():
    k = gf8()
    assert k.size() == 8
    assert k.modulus_coeffs() == [1, 1, 0, 1]


def test_fq_identities():
    k = gf8()
    assert k.zero().is_zero()
    assert not k.one().is_zero()
    assert k.one().coeffs() == [1]


def test_fq_elem():
    k = gf8()
    x = k.elem([0, 1])
    assert x.coeffs() == [0, 1]
    # Representatives of degree below 3 are already reduced modulo 1 + x + x^3.
    y = k.elem([1, 0, 1])
    assert y.coeffs() == [1, 0, 1]


def test_fq_add_sub():
    k = gf8()
    a = k.elem([1, 0, 1])
    b = k.elem([0, 1])
    # In characteristic 2, addition is coefficient-wise XOR.
    s = k.add(a, b)
    assert s.coeffs() == [1, 1, 1]


def test_fq_mul():
    k = gf8()
    a = k.elem([1, 0, 1])
    b = k.elem([0, 1])
    # (1 + x^2) * x = x + x^3 = 1 after reducing by x^3 = x + 1.
    p = k.mul(a, b)
    assert p.coeffs() == [1]


def test_fq_inv():
    k = gf8()
    a = k.elem([1, 0, 1])
    inv = k.inv(a)
    prod = k.mul(a, inv)
    assert prod.coeffs() == [1]


def test_fq_inv_zero_error():
    k = gf8()
    with pytest.raises(ZeroDivisionError):
        k.inv(k.zero())


def test_fq_pow():
    k = gf8()
    x = k.elem([0, 1])
    # The chosen modulus makes x primitive in GF(8), so x has order 2^3 - 1.
    assert k.pow(x, 7).coeffs() == [1]
    assert k.pow(x, -1).coeffs() == k.inv(x).coeffs()


def test_fq_multiplicative_order():
    k = gf8()
    x = k.elem([0, 1])
    # A primitive element of GF(8)^* has order 7.
    assert k.mul_order(x) == 7


def test_fq_primitive_elements():
    k = gf8()
    prims = k.primitive_elements()
    assert len(prims) > 0
    for p in prims:
        # Primitive means generator of the cyclic multiplicative group GF(8)^*.
        assert k.mul_order(p) == k.size() - 1


def test_fq_elements_enumeration():
    k = gf8()
    elems = k.elements(None)
    assert len(elems) == 8


def test_fq_negative_power():
    k = gf8()
    a = k.elem([0, 1])
    assert k.pow(a, -1).coeffs() == k.inv(a).coeffs()


def test_fq_element_operators():
    k = gf8()
    a = k.elem([1, 0, 1])
    b = k.elem([0, 1])
    assert (a + b).coeffs() == [1, 1, 1]
    # Same product as test_fq_mul, routed through Python operator overloads.
    assert (a * b).coeffs() == [1]
    assert (b ** 7).coeffs() == [1]
    assert (a / b).coeffs() == [1, 1, 1]


def test_fq_scalar_mixed_ops():
    k = gf8()
    b = k.elem([0, 1])
    # Integer scalars are reduced into the prime subfield F2 before arithmetic.
    assert (3 + b).coeffs() == [1, 1]
    assert (2 * b).is_zero()


def test_fq_trace():
    k = gf8()
    a = k.elem([1, 0, 1])
    # Trace sums Frobenius conjugates: a + a^2 + a^4 lands in the prime subfield.
    t = a.trace()
    assert t.coeffs() == [1]


def test_fq_norm():
    k = gf8()
    a = k.elem([1, 0, 1])
    # Norm multiplies Frobenius conjugates: a * a^2 * a^4 lands in the prime subfield.
    n = a.norm()
    assert n.coeffs() == [1]


def test_fq_permutation_actions():
    k = gf8()
    a = k.elem([0, 1])
    b = k.elem([1])

    add = k.add_perm(b)
    mul = k.mul_perm(a)
    aff = k.affine_perm(a, b)

    assert add.n() == 8
    assert mul.n() == 8
    # Composition convention is add.compose(mul)(x) = add(mul(x)) = a*x + b.
    assert add.compose(mul) == aff
    # Multiplication by zero is not a bijection, so it cannot define a permutation.
    with pytest.raises(ValueError):
        k.mul_perm(k.zero())


def test_fq_multiplicative_action_helpers():
    k = gf8()
    # GF(8)^* has 7 elements and acts cyclically on the nonzero field elements.
    assert k.multiplicative_action_subgroup_size() == 7
    orbits = k.multiplicative_action_orbits()
    assert orbits == [[0], [1, 2, 3, 4, 5, 6, 7]]
    # Nonzero points have trivial stabilizer; zero is fixed by all multiplications.
    assert k.multiplicative_action_stabilizer_size(1) == 1
    assert k.multiplicative_action_stabilizer_size(0) == 7
    assert not k.is_multiplicative_action_transitive()

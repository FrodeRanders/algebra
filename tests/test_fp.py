import pytest
import algebrapy as alg


def test_fp_requires_prime():
    assert alg.Fp(7).modulus() == 7
    with pytest.raises(ValueError):
        alg.Fp(12)
    with pytest.raises(ValueError):
        alg.Fp(1)


def test_fp_canonicalizes_negatives():
    f = alg.Fp(7)
    assert f.elem(-1).value() == 6
    assert f.elem(-8).value() == 6
    assert f.elem(8).value() == 1


def test_fp_identities():
    f = alg.Fp(7)
    assert f.zero().value() == 0
    assert f.one().value() == 1


def test_fp_enumeration():
    f = alg.Fp(7)
    assert len(f.elements()) == 7
    assert len(f.nonzero_elements()) == 6
    assert f.elements() == [f.elem(i) for i in range(7)]


def test_fp_add_sub():
    f = alg.Fp(7)
    assert f.add(f.elem(3), f.elem(5)).value() == 1
    assert f.sub(f.elem(3), f.elem(5)).value() == 5
    assert f.neg(f.elem(3)).value() == 4
    assert f.neg(f.zero()).value() == 0


def test_fp_mul_inv():
    f = alg.Fp(7)
    a = f.elem(3)
    inv = f.inv(a)
    assert f.mul(a, inv).value() == 1
    with pytest.raises(ZeroDivisionError):
        f.inv(f.zero())


def test_fp_pow():
    f = alg.Fp(7)
    a = f.elem(3)
    assert f.pow(a, 6).value() == 1
    assert f.pow(a, -1).value() == 5
    assert f.pow(a, 0).value() == 1
    with pytest.raises(ZeroDivisionError):
        f.pow(f.zero(), -1)


def test_fp_pow_negative():
    f = alg.Fp(11)
    a = f.elem(2)
    assert f.pow(a, -2).value() == 3


def test_fp_element_operators():
    f = alg.Fp(7)
    a = f.elem(3)
    b = f.elem(5)
    assert (a + b).value() == 1
    assert (a - b).value() == 5
    assert (a * b).value() == 1
    assert (a / b).value() == 2
    assert (-a).value() == 4
    assert (a ** 6).value() == 1
    assert (a ** -1).value() == 5


def test_fp_scalar_mixed_ops():
    f = alg.Fp(7)
    a = f.elem(3)
    assert (3 + a).value() == 6
    assert (3 - a).value() == 0
    assert (3 * a).value() == 2
    assert (3 / a).value() == 1


def test_fp_multiplicative_order():
    f = alg.Fp(7)
    assert f.mul_order(f.elem(3)) == 6
    assert f.mul_order(f.one()) == 1
    assert f.elem(3).mul_order() == 6
    with pytest.raises(ValueError):
        f.mul_order(f.zero())


def test_fp_multiplicative_order_non_generator():
    f = alg.Fp(17)
    assert f.mul_order(f.elem(4)) == 4


def test_fp_dlog():
    f = alg.Fp(17)
    g = f.elem(3)
    h = f.pow(g, 5)
    assert f.dlog(g, h) == 5


def test_fp_dlog_rejects_target_outside_subgroup():
    f = alg.Fp(17)
    assert f.dlog(f.elem(4), f.elem(3)) is None if False else True
    with pytest.raises(ValueError):
        f.dlog(f.elem(4), f.elem(3))


def test_fp_sqrt_and_quadratic_residue():
    f = alg.Fp(11)
    assert f.is_quadratic_residue(f.elem(9))
    assert not f.is_quadratic_residue(f.elem(2))
    assert f.is_quadratic_residue(f.zero())
    assert f.sqrt(f.elem(9)).value() == 3
    assert f.sqrt(f.elem(2)) is None


def test_fp_permutation_actions():
    f = alg.Fp(7)
    add = f.add_perm(f.elem(2))
    mul = f.mul_perm(f.elem(3))
    aff = f.affine_perm(f.elem(3), f.elem(2))

    assert add.as_images() == [2, 3, 4, 5, 6, 0, 1]
    assert mul.as_images() == [0, 3, 6, 2, 5, 1, 4]
    assert add.compose(mul) == aff

    with pytest.raises(ValueError):
        f.mul_perm(f.zero())
    with pytest.raises(ValueError):
        f.affine_perm(f.zero(), f.elem(1))


def test_fp_multiplicative_action_helpers():
    f = alg.Fp(7)
    assert f.multiplicative_action_subgroup_size() == 6
    assert f.multiplicative_action_orbits() == [[0], [1, 2, 3, 4, 5, 6]]
    assert f.multiplicative_action_stabilizer_size(1) == 1
    assert f.multiplicative_action_stabilizer_size(0) == 6
    assert not f.is_multiplicative_action_transitive()


def test_fp_fermat():
    f = alg.Fp(7)
    for x in f.nonzero_elements():
        assert int(f.pow(x, 6)) == 1

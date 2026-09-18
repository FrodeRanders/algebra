import pytest
import algebrapy as alg


def test_fp_actions_are_lazy_and_match_perms():
    f = alg.Fp(7)
    affine = f.affine_action(f.elem(3), f.elem(2))

    assert affine.size() == 7
    assert affine.kind() == "affine"
    assert affine.apply(2) == 1
    assert affine.as_perm() == f.affine_perm(f.elem(3), f.elem(2))
    assert f.add_action(f.elem(2)).kind() == "translation"
    assert f.mul_action(f.elem(3)).kind() == "multiplication"


def test_fp_action_compose_and_inverse():
    f = alg.Fp(7)
    add = f.add_action(f.elem(2))
    mul = f.mul_action(f.elem(3))
    composed = add.compose(mul)

    for x in range(7):
        assert composed.apply(x) == (3 * x + 2) % 7

    inverse = composed.inverse()
    for x in range(7):
        assert inverse.apply(composed.apply(x)) == x


def test_action_cycle_and_bounds():
    f = alg.Fp(7)
    mul = f.mul_action(f.elem(3))

    assert mul.cycle(0) == [0]
    assert mul.cycle(1) == [1, 3, 2, 6, 4, 5]
    with pytest.raises(ValueError):
        mul.apply(7)
    with pytest.raises(ValueError):
        mul.cycle(7)


def test_action_domain_mismatch_and_errors():
    f7 = alg.Fp(7)
    f11 = alg.Fp(11)

    with pytest.raises(ValueError):
        f7.add_action(f7.one()).compose(f11.add_action(f11.one()))
    with pytest.raises(ValueError):
        alg.Zn(7).add_action(alg.Zn(7).one()).compose(f7.add_action(f7.one()))
    with pytest.raises(ValueError):
        f7.mul_action(f7.zero())
    with pytest.raises(ValueError):
        f7.affine_action(f7.zero(), f7.one())


def test_zn_actions_require_units():
    z = alg.Zn(12)
    five = z.elem(5)

    assert z.mul_action(five).apply(7) == 11
    with pytest.raises(ValueError):
        z.mul_action(z.elem(6))
    with pytest.raises(ValueError):
        z.affine_action(z.elem(6), z.elem(1))

    affine = z.affine_action(five, z.elem(1))
    assert affine.apply(2) == 11
    inverse = affine.inverse()
    for x in range(12):
        assert inverse.apply(affine.apply(x)) == x


def test_fq_action_matches_perm():
    f = alg.Fq(2, [1, 1, 0, 0, 1])  # x^4 + x + 1
    a = f.elem([1, 1])
    b = f.elem([1, 0, 1])

    assert f.affine_action(a, b).as_perm(16) == f.affine_perm(a, b)
    assert f.mul_action(a).as_perm(16) == f.mul_perm(a)
    assert f.add_action(b).as_perm(16) == f.add_perm(b)


def test_fq_large_field_without_enumeration():
    # GF(2^16) with x^16 + x^12 + x^3 + x + 1.
    f = alg.Fq(2, [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1])
    assert f.size() == 65536

    shift = f.add_action(f.one())
    assert shift.apply(0) == 1
    assert shift.apply(1) == 0
    assert shift.cycle(0) == [0, 1]

    x = f.elem([0, 1])
    mul = f.mul_action(x)
    cycle = mul.cycle(1)
    assert len(cycle) == 65535
    assert 0 not in cycle
    with pytest.raises(ValueError):
        mul.as_perm(1024)

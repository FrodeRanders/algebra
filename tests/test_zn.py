import pytest
import algebrapy as alg


def test_zn_requires_modulus_at_least_two():
    assert alg.Zn(2).modulus() == 2
    with pytest.raises(ValueError):
        alg.Zn(1)


def test_zn_elements():
    z = alg.Zn(8)
    assert len(z.elements()) == 8
    assert z.zero().value() == 0
    assert z.one().value() == 1


def test_zn_canonicalizes_negatives():
    z = alg.Zn(8)
    assert z.elem(-1).value() == 7
    assert z.elem(-9).value() == 7
    assert z.elem(9).value() == 1


def test_zn_units():
    z = alg.Zn(8)
    units = [u.value() for u in z.units()]
    assert units == [1, 3, 5, 7]
    assert z.is_unit(z.elem(3))
    assert not z.is_unit(z.elem(6))


def test_zn_zero_divisors():
    z = alg.Zn(12)
    zd = [x.value() for x in z.zero_divisors()]
    assert zd == [2, 3, 4, 6, 8, 9, 10]
    assert z.is_zero_divisor(z.elem(6))
    assert not z.is_zero_divisor(z.zero())
    assert not z.is_zero_divisor(z.elem(5))


def test_zn_integral_domain():
    assert alg.Zn(7).is_integral_domain()
    assert not alg.Zn(12).is_integral_domain()


def test_zn_ideal_generation():
    z = alg.Zn(12)
    i = z.ideal(z.elem(8))
    assert i.generator().value() == 4
    assert [x.value() for x in i.elements()] == [0, 4, 8]
    assert i.size() == 3
    assert i.contains(z.elem(8))
    assert not i.contains(z.elem(6))


def test_zn_ideals_enumerated_by_divisors():
    z = alg.Zn(12)
    gens = [i.generator().value() for i in z.ideals()]
    assert gens == [1, 2, 3, 4, 6, 0]


def test_zn_ideal_prime_and_maximal():
    z = alg.Zn(12)
    i2 = z.ideal(z.elem(2))
    i3 = z.ideal(z.elem(3))
    i4 = z.ideal(z.elem(4))

    assert i2.is_prime()
    assert i2.is_maximal()
    assert i2.quotient_ring().modulus() == 2

    assert i3.is_prime()
    assert i3.is_maximal()
    assert i3.quotient_ring().modulus() == 3

    assert not i4.is_prime()
    assert not i4.is_maximal()
    assert i4.quotient_ring().modulus() == 4


def test_zero_ideal_in_prime_field_is_prime_and_maximal():
    z7 = alg.Zn(7)
    zero = z7.ideal(z7.zero())
    assert zero.is_zero()
    assert zero.is_prime()
    assert zero.is_maximal()
    assert zero.is_proper()
    assert zero.quotient_ring().modulus() == 7


def test_zn_add_sub():
    z = alg.Zn(12)
    a = z.elem(7)
    b = z.elem(11)
    assert z.add(a, b).value() == 6
    assert z.sub(a, b).value() == 8
    assert z.neg(a).value() == 5


def test_zn_mul_inv():
    z = alg.Zn(10)
    a = z.elem(3)
    inv = z.inv(a)
    assert z.mul(a, inv).value() == 1
    with pytest.raises(ValueError):
        z.inv(z.elem(2))


def test_zn_pow():
    z = alg.Zn(12)
    assert z.pow(z.elem(5), 3).value() == 5
    assert z.pow(z.elem(5), 0).value() == 1
    assert z.pow(z.elem(5), -1).value() == 5


def test_zn_negative_power_requires_unit():
    z = alg.Zn(9)
    assert z.pow(z.elem(2), -1).value() == 5
    with pytest.raises(ValueError):
        z.pow(z.elem(3), -1)


def test_zn_element_operators():
    z = alg.Zn(12)
    a = z.elem(7)
    b = z.elem(11)
    assert (a + b).value() == 6
    assert (a * b).value() == 5
    assert (-a).value() == 5
    assert (a ** 3).value() == 7
    assert a.inv().value() == 7


def test_zn_scalar_mixed_ops():
    z = alg.Zn(12)
    a = z.elem(7)
    assert (5 + a).value() == 0
    assert (5 - a).value() == 10
    assert (5 * a).value() == 11
    assert (5 / a).value() == 11


def test_zn_bool():
    z = alg.Zn(12)
    assert not z.zero().is_zero() is False
    assert z.zero().is_zero() is True
    assert bool(z.elem(5)) is True
    assert bool(z.zero()) is False


def test_zn_add_perm_translation():
    z = alg.Zn(5)
    p = z.add_perm(z.elem(2))
    assert p.as_images() == [2, 3, 4, 0, 1]
    assert p.order() == 5


def test_zn_mul_perm_exactly_for_units():
    z = alg.Zn(12)
    p = z.mul_perm(z.elem(5))
    assert p.as_images() == [0, 5, 10, 3, 8, 1, 6, 11, 4, 9, 2, 7]
    assert p.order() == 2
    with pytest.raises(ValueError):
        z.mul_perm(z.elem(6))


def test_zn_affine_perm():
    z = alg.Zn(7)
    a = z.elem(3)
    b = z.elem(2)
    mul = z.mul_perm(a)
    add = z.add_perm(b)
    aff = z.affine_perm(a, b)
    assert add.compose(mul) == aff


def test_zn_unit_group_perms():
    z = alg.Zn(10)
    perms = z.unit_group_perms()
    assert len(perms) == len(z.units())


def test_zn_unit_action_helpers():
    z = alg.Zn(12)
    assert z.unit_action_subgroup_size() == 4
    orbits = z.unit_action_orbits()
    assert orbits == [[0], [1, 5, 7, 11], [2, 10], [3, 9], [4, 8], [6]]
    assert not z.is_unit_action_transitive()
    assert z.unit_action_stabilizer_size(1) == 1
    assert z.unit_action_stabilizer_size(6) == 4


def test_zn_prime_modulus_unit_action_transitive():
    z = alg.Zn(7)
    assert z.unit_action_subgroup_size() == 6
    assert z.unit_action_orbits() == [[0], [1, 2, 3, 4, 5, 6]]
    assert not z.is_unit_action_transitive()

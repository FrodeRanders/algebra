import pytest
import algebrapy as alg


def test_perm_new():
    p = alg.Perm(3, [1, 2, 0])
    assert p.n() == 3
    assert p.as_images() == [1, 2, 0]


def test_perm_new_rejects_invalid():
    with pytest.raises(ValueError):
        alg.Perm(3, [1, 1, 2])
    with pytest.raises(ValueError):
        alg.Perm(3, [0, 1])
    with pytest.raises(ValueError):
        alg.Perm(3, [0, 1, 3])


def test_perm_identity():
    p = alg.Perm.identity(4)
    assert p.as_images() == [0, 1, 2, 3]
    assert p.order() == 1


def test_perm_cycle():
    p = alg.Perm.cycle(5, [0, 1, 2])
    assert p.as_images() == [1, 2, 0, 3, 4]


def test_perm_cycle_length_one_returns_identity():
    p = alg.Perm.cycle(4, [0])
    assert p == alg.Perm.identity(4)


def test_perm_compose():
    p = alg.Perm(4, [1, 2, 3, 0])
    q = alg.Perm(4, [0, 2, 1, 3])
    r = p.compose(q)
    assert r.as_images() == [1, 3, 2, 0]  # p(q(i))


def test_perm_mul_compose():
    p = alg.Perm(4, [1, 2, 3, 0])
    q = alg.Perm(4, [0, 2, 1, 3])
    assert (p * q).as_images() == [1, 3, 2, 0]  # p(q(i))


def test_perm_inv():
    p = alg.Perm(4, [1, 2, 3, 0])
    inv = p.inv()
    assert inv.as_images() == [3, 0, 1, 2]
    assert p.compose(inv) == alg.Perm.identity(4)


def test_perm_pow():
    p = alg.Perm.cycle(5, [0, 1, 2])
    assert p.pow(3).as_images() == [0, 1, 2, 3, 4]
    assert p.pow(-1).as_images() == [2, 0, 1, 3, 4]


def test_perm_order():
    p = alg.Perm(6, [1, 0, 2, 4, 5, 3])
    assert p.order() == 6


def test_perm_cycles():
    p = alg.Perm(6, [1, 0, 2, 4, 5, 3])
    assert p.cycles() == [[0, 1], [3, 4, 5]]


def test_perm_cycle_type():
    p = alg.Perm(6, [1, 0, 2, 4, 5, 3])
    assert p.cycle_type() == [2, 3]


def test_perm_cycle_notation():
    p = alg.Perm(6, [1, 0, 2, 4, 5, 3])
    assert p.cycle_notation() == "(0 1)(3 4 5)"
    assert alg.Perm.identity(4).cycle_notation() == "()"


def test_perm_conjugate():
    p = alg.Perm.cycle(5, [0, 1, 2])
    g = alg.Perm.cycle(5, [0, 4, 3, 2, 1])
    conj = p.conjugate_by(g)
    assert conj.cycle_type() == p.cycle_type()


def test_sn_identity():
    s4 = alg.Sn(4)
    assert s4.identity() == alg.Perm.identity(4)


def test_sn_elements():
    s3 = alg.Sn(3)
    assert len(s3.elements()) == 6
    s4 = alg.Sn(4)
    assert len(s4.elements()) == 24


def test_sn_generated_subgroup():
    s4 = alg.Sn(4)
    t = alg.Perm.cycle(4, [0, 1])
    c = alg.Perm.cycle(4, [0, 1, 2, 3])
    g = s4.generated([t, c])
    assert g.order() == 24
    assert not g.is_abelian()


def test_perm_subgroup_order_and_elements():
    s3 = alg.Sn(3)
    sigma = alg.Perm.cycle(3, [0, 1, 2])
    tau = alg.Perm.cycle(3, [0, 1])
    h = s3.generated([sigma, tau])
    assert h.order() == 6
    assert len(h.elements()) == 6


def test_perm_subgroup_is_abelian():
    s3 = alg.Sn(3)
    sigma = alg.Perm.cycle(3, [0, 1, 2])
    tau = alg.Perm.cycle(3, [0, 1])
    h = s3.generated([sigma, tau])
    assert not h.is_abelian()

    g = s3.generated([sigma])
    assert g.is_abelian()


def test_perm_subgroup_contains():
    s4 = alg.Sn(4)
    g = s4.generated([alg.Perm.cycle(4, [0, 1, 2, 3])])
    assert g.contains(alg.Perm.identity(4))
    assert not g.contains(alg.Perm.cycle(4, [0, 1]))


def test_orbit():
    s5 = alg.Sn(5)
    g = alg.Perm.cycle(5, [1, 3, 4])
    orbit = s5.orbit(1, [g])
    assert orbit == [1, 3, 4]


def test_orbits():
    s6 = alg.Sn(6)
    a = alg.Perm.cycle(6, [0, 1, 2])
    b = alg.Perm.cycle(6, [3, 4])
    orbits = s6.orbits([a, b])
    orbits.sort()
    assert orbits == [[0, 1, 2], [3, 4], [5]]


def test_subgroup_transitivity():
    s4 = alg.Sn(4)
    g = alg.Perm.cycle(4, [0, 1, 2, 3])
    assert s4.subgroup_size([g]) == 4
    assert s4.is_transitive([g])


def test_stabilizer_orbit_stabilizer_theorem():
    s4 = alg.Sn(4)
    a = alg.Perm.cycle(4, [0, 1, 2, 3])
    b = alg.Perm.cycle(4, [1, 3])
    gens = [a, b]

    group = s4.generated(gens)
    subgroup_size = group.order()
    orbit_size = len(s4.orbit(0, gens))
    stabilizer_size = s4.stabilizer_size(0, gens)
    assert subgroup_size == orbit_size * stabilizer_size


def test_conjugacy_class():
    s4 = alg.Sn(4)
    t = alg.Perm.cycle(4, [0, 1])
    size = s4.conjugacy_class_size(t, 24)
    assert size == 6

    c3 = alg.Perm.cycle(4, [0, 1, 2])
    size3 = s4.conjugacy_class_size(c3, 24)
    assert size3 == 8


def test_conjugation_preserves_cycle_type():
    p = alg.Perm.cycle(5, [0, 1, 2])
    g = alg.Perm.cycle(5, [0, 4, 3, 2, 1])
    assert p.conjugate_by(g).cycle_type() == p.cycle_type()


def test_subgroup_conjugation_and_intersection():
    s4 = alg.Sn(4)
    h = s4.generated([alg.Perm.cycle(4, [0, 1])])
    g = alg.Perm.cycle(4, [0, 2])
    conj = h.conjugate_by(g)
    expected = alg.Perm.cycle(4, [1, 2])
    assert conj.contains(expected)
    assert h.intersection(conj).order() == 1


def test_normality():
    s4 = alg.Sn(4)
    v4 = s4.generated([
        alg.Perm(4, [1, 0, 3, 2]),
        alg.Perm(4, [2, 3, 0, 1]),
    ])
    s4_full = s4.generated([
        alg.Perm.cycle(4, [0, 1]),
        alg.Perm.cycle(4, [0, 1, 2, 3]),
    ])
    assert v4.is_normal_in(s4_full)


def test_sylow_subgroups():
    s3 = alg.Sn(3)
    sigma = alg.Perm.cycle(3, [0, 1, 2])
    tau = alg.Perm.cycle(3, [0, 1])
    g = s3.generated([sigma, tau])

    sylow2 = g.sylow_p_subgroups(2)
    sylow3 = g.sylow_p_subgroups(3)

    assert len(sylow2) == 3
    assert all(h.order() == 2 for h in sylow2)
    assert len(sylow3) == 1
    assert sylow3[0].order() == 3
    assert sylow3[0].is_normal_in(g)


def test_sylow_subgroups_s4():
    s4 = alg.Sn(4)
    g = s4.generated([
        alg.Perm.cycle(4, [0, 1]),
        alg.Perm.cycle(4, [0, 1, 2, 3]),
    ])

    sylow2 = g.sylow_p_subgroups(2)
    sylow3 = g.sylow_p_subgroups(3)

    assert len(sylow2) == 3
    assert all(h.order() == 8 for h in sylow2)
    assert len(sylow3) == 4
    assert all(h.order() == 3 for h in sylow3)
    assert all(h.is_sylow_p_subgroup(2, g) for h in sylow2)


def test_dihedral_8():
    s4 = alg.Sn(4)
    r = alg.Perm.cycle(4, [0, 1, 2, 3])
    f = alg.Perm(4, [0, 3, 2, 1])
    group = s4.generated([r, f])

    assert group.contains(alg.Perm.identity(4))
    assert not group.is_abelian()
    assert group.p_part_order(2) == 8
    assert group.is_p_group(2)
    assert not group.is_p_group(3)

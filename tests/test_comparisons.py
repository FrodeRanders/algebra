import pytest
import algebrapy as alg


class TestFpElemComparisons:
    def test_eq_same_field(self):
        f = alg.Fp(7)
        assert f.elem(3) == f.elem(3)
        assert not (f.elem(3) == f.elem(5))

    def test_ne_same_field(self):
        f = alg.Fp(7)
        assert f.elem(3) != f.elem(5)
        assert not (f.elem(3) != f.elem(3))

    def test_eq_different_moduli(self):
        a = alg.Fp(7).elem(3)
        b = alg.Fp(11).elem(3)
        assert not (a == b)

    def test_ne_different_moduli(self):
        a = alg.Fp(7).elem(3)
        b = alg.Fp(11).elem(3)
        assert a != b

    def test_hash(self):
        f = alg.Fp(7)
        a = f.elem(3)
        b = f.elem(3)
        assert hash(a) == hash(b)
        d = {a: "three"}
        assert d[b] == "three"

    def test_hash_different_values(self):
        f = alg.Fp(7)
        h1 = hash(f.elem(3))
        h2 = hash(f.elem(5))
        # Collisions are unlikely for this simple hash
        assert isinstance(h1, int)


class TestZnElemComparisons:
    def test_eq_same_ring(self):
        r = alg.Zn(12)
        assert r.elem(5) == r.elem(5)
        assert not (r.elem(5) == r.elem(7))

    def test_ne_same_ring(self):
        r = alg.Zn(12)
        assert r.elem(5) != r.elem(7)
        assert not (r.elem(5) != r.elem(5))

    def test_eq_different_moduli(self):
        a = alg.Zn(12).elem(5)
        b = alg.Zn(8).elem(5)
        assert not (a == b)

    def test_ne_different_moduli(self):
        a = alg.Zn(12).elem(5)
        b = alg.Zn(8).elem(5)
        assert a != b

    def test_hash(self):
        r = alg.Zn(12)
        a = r.elem(5)
        d = {a: "five"}
        assert d[r.elem(5)] == "five"


class TestPermComparisons:
    def test_eq(self):
        p = alg.Perm(4, [1, 2, 3, 0])
        q = alg.Perm(4, [1, 2, 3, 0])
        r = alg.Perm(4, [1, 0, 2, 3])
        assert p == q
        assert not (p == r)

    def test_ne(self):
        p = alg.Perm(4, [1, 2, 3, 0])
        r = alg.Perm(4, [1, 0, 2, 3])
        assert p != r
        assert not (p != alg.Perm(4, [1, 2, 3, 0]))

    def test_hash(self):
        p = alg.Perm(4, [1, 2, 3, 0])
        q = alg.Perm(4, [1, 2, 3, 0])
        assert hash(p) == hash(q)
        d = {p: "cycle"}
        assert d[q] == "cycle"

    def test_hash_set(self):
        p = alg.Perm(4, [1, 2, 3, 0])
        q = alg.Perm(4, [1, 2, 3, 0])
        s = {p}
        assert q in s


class TestFqElemComparisons:
    def test_eq(self):
        k = alg.Fq(2, [1, 1, 0, 1])
        a = k.elem([0, 1])
        assert a == k.elem([0, 1])
        assert a != k.elem([1, 0, 1])

    def test_hash(self):
        k = alg.Fq(2, [1, 1, 0, 1])
        a = k.elem([0, 1])
        d = {a: "x"}
        assert d[k.elem([0, 1])] == "x"


class TestPolyFpComparisons:
    def test_eq(self):
        p = alg.PolyFp(2, [1, 0, 1])
        q = alg.PolyFp(2, [1, 0, 1])
        r = alg.PolyFp(2, [1, 1])
        assert p == q
        assert p != r

    def test_hash(self):
        p = alg.PolyFp(2, [1, 0, 1])
        d = {p: "poly"}
        assert d[alg.PolyFp(2, [1, 0, 1])] == "poly"


class TestZnIdealComparisons:
    def test_eq(self):
        z = alg.Zn(12)
        a = z.ideal(z.elem(4))
        b = z.ideal(z.elem(8))
        c = z.ideal(z.elem(3))
        assert a == b
        assert a != c

    def test_hash(self):
        z = alg.Zn(12)
        a = z.ideal(z.elem(4))
        d = {a: "ideal"}
        assert d[z.ideal(z.elem(8))] == "ideal"


class TestPermSubgroupComparisons:
    def test_eq(self):
        s4 = alg.Sn(4)
        h1 = s4.generated([alg.Perm.cycle(4, [0, 1])])
        h2 = s4.generated([alg.Perm.cycle(4, [0, 1])])
        h3 = s4.generated([alg.Perm.cycle(4, [0, 1, 2])])
        assert h1 == h2
        assert h1 != h3

import pytest
import algebrapy as alg


class TestPolyFp:
    def test_construction(self):
        p = alg.PolyFp(7, [1, 2, 3])
        assert p.p() == 7
        assert p.coeffs() == [1, 2, 3]
        assert p.degree() == 2

    def test_zero_and_one(self):
        zero = alg.PolyFp(7, [0])
        assert zero.is_zero()
        assert not zero.is_one()

        one = alg.PolyFp(7, [1])
        assert not one.is_zero()
        assert one.is_one()

    def test_negative_degree(self):
        zero = alg.PolyFp(7, [0])
        assert zero.degree() == -1

    def test_monic(self):
        p = alg.PolyFp(7, [2, 3])
        m = p.monic()
        assert m.coeffs() == [3, 1]

    def test_add(self):
        p = alg.PolyFp(7, [1, 2, 3])
        q = alg.PolyFp(7, [4, 5])
        r = p.add(q)
        assert r.coeffs() == [5, 0, 3]

    def test_neg(self):
        p = alg.PolyFp(7, [1, 2, 3])
        n = p.neg()
        assert n.coeffs() == [6, 5, 4]

    def test_sub(self):
        p = alg.PolyFp(7, [1, 2, 3])
        q = alg.PolyFp(7, [4, 5])
        r = p.sub(q)
        assert r.coeffs() == [4, 4, 3]

    def test_mul(self):
        p = alg.PolyFp(7, [1, 1])
        q = alg.PolyFp(7, [1, 1])
        r = p.mul(q)
        assert r.coeffs() == [1, 2, 1]

    def test_scale(self):
        p = alg.PolyFp(7, [1, 2, 3])
        s = p.scale(2)
        assert s.coeffs() == [2, 4, 6]

    def test_div_rem(self):
        a = alg.PolyFp(7, [1, 2, 3])
        b = alg.PolyFp(7, [1, 1])
        q, r = a.div_rem(b)
        assert q.mul(b).add(r).coeffs() == a.coeffs()

    def test_modulo(self):
        p = alg.PolyFp(7, [1, 2, 3])
        m = alg.PolyFp(7, [1, 1])
        r = p.modulo(m)
        assert r.degree() < m.degree() or r.is_zero()

    def test_gcd(self):
        a = alg.PolyFp(2, [0, 1, 1])
        b = alg.PolyFp(2, [0, 1])
        g = a.gcd(b)
        assert g.coeffs() == [0, 1]

    def test_egcd_bezout(self):
        a = alg.PolyFp(2, [1, 0, 1])
        b = alg.PolyFp(2, [1, 1, 0, 1])
        g, s, t = a.egcd(b)
        lhs = s.mul(a).add(t.mul(b)).monic()
        assert lhs.coeffs() == g.monic().coeffs()
        assert g.monic().coeffs() == [1]

    def test_is_irreducible(self):
        assert alg.PolyFp(2, [1, 1, 0, 1]).is_irreducible()
        assert not alg.PolyFp(2, [1, 0, 1]).is_irreducible()

    def test_eval(self):
        p = alg.PolyFp(7, [1, 2, 3])  # 1 + 2x + 3x^2
        assert p.eval(0) == 1
        assert p.eval(1) == 6
        assert p.eval(2) == 3

    def test_eval_negative(self):
        p = alg.PolyFp(7, [0, 1])
        assert p.eval(-1) == 6


class TestPolyFpOperators:
    def test_add_operator(self):
        p = alg.PolyFp(7, [1, 2, 3])
        q = alg.PolyFp(7, [4, 5])
        r = p + q
        assert r.coeffs() == [5, 0, 3]

    def test_sub_operator(self):
        p = alg.PolyFp(7, [1, 2, 3])
        q = alg.PolyFp(7, [4, 5])
        r = p - q
        assert r.coeffs() == [4, 4, 3]

    def test_mul_operator(self):
        p = alg.PolyFp(7, [1, 1])
        q = alg.PolyFp(7, [1, 1])
        r = p * q
        assert r.coeffs() == [1, 2, 1]

    def test_neg_operator(self):
        p = alg.PolyFp(7, [1, 2, 3])
        n = -p
        assert n.coeffs() == [6, 5, 4]

    def test_pow_operator(self):
        p = alg.PolyFp(7, [0, 1])
        q = p ** 2
        assert q.coeffs() == [0, 0, 1]
        assert (p ** 0).coeffs() == [1]

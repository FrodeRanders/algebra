import pytest
import algebrapy as alg


class TestCn:
    def test_construction(self):
        c = alg.Cn(12)
        assert c.order() == 12
        assert len(c) == 12
        assert c.identity() == 0

    def test_requires_valid_order(self):
        with pytest.raises(ValueError):
            alg.Cn(0)

    def test_trivial_group(self):
        c = alg.Cn(1)
        assert c.elements() == [0]
        assert c.generators() == []
        assert c.element_order(0) == 1

    def test_elements(self):
        c = alg.Cn(12)
        assert c.elements() == list(range(12))

    def test_generators_phi12(self):
        c = alg.Cn(12)
        assert c.generators() == [1, 5, 7, 11]
        assert c.is_generator(5)
        assert not c.is_generator(4)
        assert not c.is_generator(12)

    def test_element_order(self):
        c = alg.Cn(12)
        assert c.element_order(0) == 1
        assert c.element_order(4) == 3
        assert c.element_order(6) == 2
        assert c.element_order(5) == 12

    def test_element_order_out_of_range(self):
        c = alg.Cn(12)
        with pytest.raises(ValueError):
            c.element_order(12)

    def test_op(self):
        c = alg.Cn(12)
        assert c.op(7, 8) == 3

    def test_op_out_of_range(self):
        c = alg.Cn(12)
        with pytest.raises(ValueError):
            c.op(12, 1)

    def test_inv(self):
        c = alg.Cn(12)
        assert c.inv(0) == 0
        assert c.inv(5) == 7
        assert c.inv(11) == 1

    def test_pow(self):
        c = alg.Cn(12)
        assert c.pow(5, 3) == 3
        assert c.pow(5, -1) == 7
        assert c.pow(0, 5) == 0

    def test_subgroup(self):
        c = alg.Cn(12)
        gens, elems = c.subgroup([8])
        assert gens == [4]
        assert elems == [0, 4, 8]

    def test_subgroup_empty(self):
        c = alg.Cn(12)
        gens, elems = c.subgroup([])
        assert gens == [0]
        assert elems == [0]

    def test_subgroup_whole_group(self):
        c = alg.Cn(12)
        gens, elems = c.subgroup([5])
        assert gens == [1]
        assert len(elems) == 12

    def test_subgroups(self):
        c = alg.Cn(12)
        subs = c.subgroups()
        assert len(subs) == 6
        sizes = [len(elems) for _, elems in subs]
        assert sizes == [12, 6, 4, 3, 2, 1]

    def test_repr(self):
        c = alg.Cn(12)
        assert repr(c) == "C12"

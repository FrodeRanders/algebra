import pytest
import algebrapy as alg


def hamming():
    return alg.BinaryBchCode(3, [1, 1, 0, 1], 3)


def bch15():
    return alg.BinaryBchCode(4, [1, 1, 0, 0, 1], 5)


def rs8():
    return alg.ReedSolomonCode(2, [1, 1, 0, 1], 3)


class TestBinaryBch:
    def test_hamming_construction(self):
        code = hamming()
        assert code.length() == 7
        assert code.dimension() == 4
        assert code.generator_poly().coeffs() == [1, 1, 0, 1]
        assert code.correction_capacity() == 1

    def test_bch15_construction(self):
        code = bch15()
        assert code.length() == 15
        assert code.generator_degree() == 8
        assert code.dimension() == 7
        assert code.correction_capacity() == 2
        assert len(code.generator_poly().coeffs()) == 9

    def test_encode_produces_divisible_codeword(self):
        code = hamming()
        msg = [1, 0, 1, 1]
        word = code.encode(msg)
        assert len(word) == 7
        assert code.is_codeword(word)
        assert code.parity_check(word)
        assert code.extract_message(word) == msg

    def test_systematic_round_trip(self):
        code = hamming()
        msg = [1, 0, 1, 1]
        word = code.encode_systematic(msg)
        assert code.is_codeword(word)
        assert code.extract_systematic_message(word) == msg

    def test_parity_check_detects_error(self):
        code = hamming()
        word = code.encode([1, 1, 0, 1])
        word[2] ^= 1
        assert not code.is_codeword(word)
        assert not code.parity_check(word)
        assert any(not (len(s) == 0) for s in code.syndromes(word))

    def test_hamming_decode_single_error(self):
        code = hamming()
        msg = [1, 0, 1, 1]
        sent = code.encode(msg)
        received = sent.copy()
        received[4] ^= 1
        corrected = code.decode(received)
        assert corrected == sent
        assert code.parity_check(corrected)
        assert code.decode_message(received) == msg

    def test_systematic_decode_after_error(self):
        code = hamming()
        msg = [1, 0, 1, 1]
        received = code.encode_systematic(msg)
        received[5] ^= 1
        assert code.decode_systematic_message(received) == msg

    def test_hamming_decode_all_single_errors(self):
        code = hamming()
        for mask in range(1 << code.dimension()):
            msg = _bits_from_mask(mask, code.dimension())
            codeword = code.encode_systematic(msg)
            for pos in range(code.length()):
                received = codeword.copy()
                received[pos] ^= 1
                assert code.decode(received) == codeword
                assert code.decode_systematic_message(received) == msg

    def test_bch15_decode_two_errors(self):
        code = bch15()
        msg = [1, 0, 1, 1, 0, 1, 1]
        sent = code.encode(msg)
        received = sent.copy()
        received[2] ^= 1
        received[11] ^= 1
        corrected = code.decode(received)
        assert corrected == sent
        assert code.extract_message(sent) == msg

    def test_decode_returns_codeword_on_three_errors(self):
        code = bch15()
        sent = code.encode([1, 0, 1, 1, 0, 1, 1])
        received = sent.copy()
        received[1] ^= 1
        received[4] ^= 1
        received[7] ^= 1
        decoded = code.decode(received)
        assert code.is_codeword(decoded)

    def test_extract_message_rejects_non_codeword(self):
        code = hamming()
        with pytest.raises(ValueError):
            code.extract_message([1, 0, 0, 0, 0, 0, 0])

    def test_shortened_parameters(self):
        code = bch15()
        n_short, k_short = code.shortened_parameters(3)
        assert n_short == 12
        assert k_short == 4

    def test_shortened_round_trip(self):
        code = bch15()
        msg = [1, 0, 1, 1]
        word = code.encode_shortened_systematic(msg, 3)
        assert len(word) == 12
        assert code.extract_shortened_systematic_message(word, 3) == msg

    def test_shortened_decode(self):
        code = bch15()
        msg = [1, 0, 1, 1]
        received = code.encode_shortened_systematic(msg, 3)
        received[2] ^= 1
        received[8] ^= 1
        assert code.decode_shortened_systematic_message(received, 3) == msg

    def test_generator_parity_check_orthogonality(self):
        code = bch15()
        g = code.generator_matrix()
        h = code.parity_check_matrix()
        h_rows = len(h)
        for row in g:
            result = [0] * h_rows
            for j, bit in enumerate(row):
                if bit:
                    for k in range(h_rows):
                        result[k] ^= h[k][j]
            assert all(b == 0 for b in result)

    def test_codeword_enumeration(self):
        code = hamming()
        codewords = code.codewords()
        assert len(codewords) == 16
        assert code.minimum_distance() == 3

    def test_weight_distribution(self):
        code = hamming()
        dist = code.weight_distribution()
        assert dist[0] == 1
        assert dist[3] == 7
        assert dist[4] == 7
        assert dist[7] == 1
        assert sum(dist) == 16


class TestReedSolomon:
    def test_parameters(self):
        rs = rs8()
        assert rs.length() == 7
        assert rs.dimension() == 3
        assert rs.designed_distance() == 5
        assert rs.correction_capacity() == 2
        assert len(rs.generator_poly()) == 5

    def test_systematic_round_trip(self):
        rs = rs8()
        fq = rs.field()
        msg = [fq.elem([1]), fq.elem([0, 1]), fq.elem([1, 1])]
        codeword = rs.encode_systematic(msg)
        assert rs.is_codeword(codeword)
        recovered = rs.extract_systematic_message(codeword)
        assert [e.coeffs() for e in recovered] == [e.coeffs() for e in msg]

    def test_decodes_two_symbol_errors(self):
        rs = rs8()
        fq = rs.field()
        msg = [fq.elem([1]), fq.elem([0, 1]), fq.elem([1, 1])]
        received = rs.encode_systematic(msg)
        received[1] = fq.add(received[1], fq.elem([1]))
        received[5] = fq.add(received[5], fq.elem([1, 1]))
        decoded = rs.decode(received)
        recovered = rs.extract_systematic_message(decoded)
        assert [e.coeffs() for e in recovered] == [e.coeffs() for e in msg]

    def test_rejects_too_many_errors(self):
        rs = rs8()
        fq = rs.field()
        msg = [fq.elem([1]), fq.elem([0, 1]), fq.elem([1, 1])]
        received = rs.encode_systematic(msg)
        received[0] = fq.add(received[0], fq.elem([1]))
        received[1] = fq.add(received[1], fq.elem([1, 1]))
        received[2] = fq.add(received[2], fq.elem([0, 1]))
        with pytest.raises(ValueError):
            rs.decode(received)

    def test_syndromes_vanish_on_codeword(self):
        rs = rs8()
        fq = rs.field()
        msg = [fq.elem([1]), fq.zero(), fq.one()]
        codeword = rs.encode_systematic(msg)
        assert all(s.is_zero() for s in rs.syndromes(codeword))


def _bits_from_mask(mask, length):
    return [(mask >> i) & 1 for i in range(length)]

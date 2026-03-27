import algebrapy as alg


def print_matrix(name, rows):
    print(name)
    for row in rows:
        print(" ", row)


def main():
    print("== Hamming / BCH(7,4,3) ==")
    hamming = alg.BinaryBchCode(3, [1, 1, 0, 1], 3)
    print(hamming)
    print("length:", hamming.length())
    print("dimension:", hamming.dimension())
    print("t:", hamming.correction_capacity())
    print("generator polynomial:", hamming.generator_poly())

    msg = [1, 0, 1, 1]
    codeword = hamming.encode(msg)
    systematic = hamming.encode_systematic(msg)
    print("message:", msg)
    print("cyclic codeword:", codeword)
    print("systematic codeword:", systematic)
    print("extract message:", hamming.extract_message(codeword))
    print("extract systematic message:", hamming.extract_systematic_message(systematic))

    received = systematic.copy()
    received[5] ^= 1
    print("received with one error:", received)
    print("syndromes:", hamming.syndromes(received))
    print("decoded codeword:", hamming.decode(received))
    print("decoded message:", hamming.decode_systematic_message(received))

    print("parity check:", hamming.parity_check(received))
    print("is decoded word a codeword?", hamming.is_codeword(hamming.decode(received)))

    print_matrix("G:", hamming.generator_matrix())
    print_matrix("H:", hamming.parity_check_matrix())
    print("minimum distance:", hamming.minimum_distance())
    print("weight distribution:", hamming.weight_distribution())
    print("number of codewords:", len(hamming.codewords()))

    print()
    print("== BCH(15,7,5) ==")
    bch = alg.BinaryBchCode(4, [1, 1, 0, 0, 1], 5)
    print(bch)
    print("length:", bch.length())
    print("dimension:", bch.dimension())
    print("t:", bch.correction_capacity())
    print("generator polynomial:", bch.generator_poly())

    msg2 = [1, 0, 1, 1, 0, 1, 1]
    codeword2 = bch.encode_systematic(msg2)
    print("message:", msg2)
    print("systematic codeword:", codeword2)

    received2 = codeword2.copy()
    received2[2] ^= 1
    received2[11] ^= 1
    print("received with two errors:", received2)
    print("syndromes:", bch.syndromes(received2))
    print("decoded codeword:", bch.decode(received2))
    print("decoded message:", bch.decode_systematic_message(received2))
    print("minimum distance:", bch.minimum_distance())
    print("weight distribution:", bch.weight_distribution())
    print("number of codewords:", len(bch.codewords()))

    print()
    print("== Shortened BCH(12,4,>=5) from BCH(15,7,5) ==")
    n_short, k_short = bch.shortened_parameters(3)
    print("shortened parameters:", (n_short, k_short))
    msg3 = [1, 0, 1, 1]
    short_codeword = bch.encode_shortened_systematic(msg3, 3)
    print("shortened message:", msg3)
    print("shortened codeword:", short_codeword)
    print(
        "shortened extracted message:",
        bch.extract_shortened_systematic_message(short_codeword, 3),
    )
    received3 = short_codeword.copy()
    received3[2] ^= 1
    received3[8] ^= 1
    print("shortened received with two errors:", received3)
    print(
        "shortened decoded message:",
        bch.decode_shortened_systematic_message(received3, 3),
    )


if __name__ == "__main__":
    main()

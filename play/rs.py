import algebrapy as alg


def main():
    rs = alg.ReedSolomonCode(2, [1, 1, 0, 1], 3)  # GF(2^3), RS(7,3,5)
    print(rs)
    print("n =", rs.length(), "k =", rs.dimension(), "d =", rs.designed_distance())
    print("t =", rs.correction_capacity())

    F = rs.field()
    one = F.elem([1])
    x = F.elem([0, 1])
    one_plus_x = F.elem([1, 1])

    print("generator polynomial:")
    for coeff in rs.generator_poly():
        print(" ", coeff)

    msg = [one, x, one_plus_x]
    codeword = rs.encode_systematic(msg)
    print("message:", msg)
    print("codeword:", codeword)
    print("is codeword?", rs.is_codeword(codeword))
    print("syndromes:", rs.syndromes(codeword))
    print("recovered message:", rs.extract_systematic_message(codeword))

    received = codeword.copy()
    received[1] = received[1] + one
    received[5] = received[5] + one_plus_x
    print("received with two symbol errors:", received)
    print("error syndromes:", rs.syndromes(received))
    print("decoded codeword:", rs.decode(received))
    print("decoded message:", rs.decode_message(received))


if __name__ == "__main__":
    main()

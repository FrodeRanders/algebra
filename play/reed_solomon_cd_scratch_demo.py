import algebrapy as alg


def interleave(codewords):
    n = len(codewords[0])
    out = []
    for col in range(n):
        for row in range(len(codewords)):
            out.append(codewords[row][col])
    return out


def deinterleave(stream, rows, cols):
    out = [[None for _ in range(cols)] for _ in range(rows)]
    idx = 0
    for col in range(cols):
        for row in range(rows):
            out[row][col] = stream[idx]
            idx += 1
    return out


def add_burst_error(stream, start, length, error_symbol):
    damaged = stream.copy()
    for i in range(start, min(start + length, len(damaged))):
        damaged[i] = damaged[i] + error_symbol
    return damaged


def count_symbol_differences(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def main():
    print("Toy CD-style Reed-Solomon scratch demo")
    print("This is not the exact CD CIRC format.")
    print("It demonstrates the same principle: interleaving spreads a burst error.")
    print()

    rs = alg.ReedSolomonCode(2, [1, 0, 0, 1, 1], 11)  # GF(16), RS(15,11,5)
    F = rs.field()
    t = rs.correction_capacity()

    print(rs)
    print("n =", rs.length(), "k =", rs.dimension(), "d =", rs.designed_distance(), "t =", t)
    print()

    messages = [
        [F.elem([1]), F.elem([0, 1]), F.elem([1, 1]), F.elem([0, 0, 1]), F.elem([1, 0, 1]),
         F.elem([1, 1, 1]), F.elem([0, 1, 1]), F.elem([1, 1, 0]), F.elem([0, 0, 0, 1]),
         F.elem([1, 0, 0, 1]), F.elem([0, 1, 0, 1])],
        [F.elem([0]), F.elem([1]), F.elem([0, 1]), F.elem([1, 1]), F.elem([0, 0, 1]),
         F.elem([1, 0, 1]), F.elem([1, 1, 1]), F.elem([0, 1, 1]), F.elem([1, 1, 0]),
         F.elem([0, 0, 0, 1]), F.elem([1, 0, 0, 1])],
        [F.elem([1, 1]), F.elem([0, 0, 1]), F.elem([1, 0, 1]), F.elem([1, 1, 1]),
         F.elem([0, 1, 1]), F.elem([1, 1, 0]), F.elem([0, 0, 0, 1]), F.elem([1, 0, 0, 1]),
         F.elem([0, 1, 0, 1]), F.elem([1, 1, 0, 1]), F.elem([0, 0, 1, 1])],
        [F.elem([0, 1]), F.elem([1, 1]), F.elem([0, 0, 1]), F.elem([1, 0, 1]), F.elem([1, 1, 1]),
         F.elem([0, 1, 1]), F.elem([1, 1, 0]), F.elem([0, 0, 0, 1]), F.elem([1, 0, 0, 1]),
         F.elem([0, 1, 0, 1]), F.elem([1, 1, 0, 1])],
    ]

    codewords = [rs.encode_systematic(msg) for msg in messages]
    error_symbol = F.elem([1])
    burst_start = 8
    burst_len = 8

    print("Scenario A: no interleaving")
    direct = codewords[0]
    direct_damaged = add_burst_error(direct, burst_start, burst_len, error_symbol)
    direct_errors = count_symbol_differences(direct, direct_damaged)
    print("symbol errors in one codeword =", direct_errors)
    try:
        rs.decode(direct_damaged)
        print("decoder succeeded")
    except Exception as exc:
        print("decoder failed:", exc)
    print()

    print("Scenario B: interleave 4 codewords, then apply the same burst to the stream")
    stream = interleave(codewords)
    damaged_stream = add_burst_error(stream, burst_start, burst_len, error_symbol)
    damaged_codewords = deinterleave(damaged_stream, len(codewords), len(codewords[0]))

    per_codeword_errors = [
        count_symbol_differences(original, damaged)
        for original, damaged in zip(codewords, damaged_codewords)
    ]
    print("symbol errors per codeword after deinterleaving =", per_codeword_errors)

    recovered_messages = []
    for i, damaged in enumerate(damaged_codewords):
        decoded = rs.decode(damaged)
        recovered = rs.extract_systematic_message(decoded)
        recovered_messages.append(recovered)
        print(f"codeword {i}: decoded successfully =", recovered == messages[i])

    print()
    print("Takeaway:")
    print("without interleaving, the burst creates too many errors in one codeword")
    print("with interleaving, the same burst is spread out and each codeword stays within t")


if __name__ == "__main__":
    main()

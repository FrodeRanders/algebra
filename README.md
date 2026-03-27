# Algebra for Python

This is algebrapy - a Python package for abstract algebra implemented in Rust (using PyO3).

## What it provides

1. Finite Fields
- Fp / FpElem - Prime fields GF(p) where p is prime
- Fq / FqElem - Extension fields GF(p^k) constructed as polynomial quotient rings Fpx/(f) where f is irreducible
2. Permutation Groups
- Perm - Permutations (bijections) with composition, inversion, exponentiation
- Sn - Symmetric groups S_n with element enumeration and subgroup generation
3. Arithmetic Utilities
- Extended GCD algorithm
- Prime number checking
- Modular inverse
4. Coding Theory
- BinaryBchCode - primitive narrow-sense binary BCH generator construction, encoding, and syndrome checks
- ReedSolomonCode - primitive narrow-sense Reed-Solomon codes over GF(p^m)

## Example usage

```python
import algebrapy as alg

# GF(7) - prime field
a = alg.Fp(7).elem(3)    # 3 (mod 7)

# GF(2^3) - extension field with irreducible polynomial 1 + x + x^3
K = alg.Fq(2, [1,1,0,1])
a = K.elem([1,0,1])      # 1 + x^2
```

The implementation supports addition, subtraction, multiplication, division, exponentiation (including negative exponents), multiplicative order, and other algebraic operations. Built with maturin for easy Python integration.

## BCH playground

```python
import algebrapy as alg

# Primitive binary BCH / Hamming code of length 7 over GF(2^3)
code = alg.BinaryBchCode(3, [1, 1, 0, 1], 3)   # x^3 + x + 1
print(code)
print("n =", code.length(), "k =", code.dimension())
print("g(x) =", code.generator_poly())

msg = [1, 0, 1, 1]
c = code.encode(msg)
print("codeword =", c)
print("is codeword?", code.is_codeword(c))
print("syndromes =", code.syndromes(c))
print("extracted message =", code.extract_message(c))

cs = code.encode_systematic(msg)
print("systematic codeword =", cs)
print("systematic message =", code.extract_systematic_message(cs))
print("G =", code.generator_matrix())
print("H =", code.parity_check_matrix())

r = c.copy()
r[2] ^= 1
print("single-bit error:", r)
print("parity check =", code.parity_check(r))
print("syndromes =", code.syndromes(r))
print("t =", code.correction_capacity())
print("decoded =", code.decode(r))
print("decoded message =", code.decode_message(r))

# Shorten BCH(15,7,5) by fixing three message coordinates to zero
big = alg.BinaryBchCode(4, [1, 1, 0, 0, 1], 5)
print("shortened parameters =", big.shortened_parameters(3))
smsg = [1, 0, 1, 1]
sc = big.encode_shortened_systematic(smsg, 3)
print("shortened codeword =", sc)
print("shortened decoded message =", big.decode_shortened_systematic_message(sc, 3))
```

There is also a fuller runnable example in [`play/bch.py`](/Users/froran/Projects/gautelis/algebra/play/bch.py) that exercises:
- cyclic and systematic encoding
- shortened systematic encoding
- syndrome computation and decoding
- message extraction
- generator and parity-check matrices
- codeword enumeration, minimum distance, and weight distribution

Current decoder scope: `BinaryBchCode.decode()` uses a binary BCH syndrome decoder based on Berlekamp-Massey and Chien search. It is intended to correct up to `t = floor((d - 1) / 2)` errors. As with standard bounded-distance BCH decoding, words with more than `t` errors can still miscorrect to a different valid codeword, so decoding beyond `t` should not be treated as reliable detection.

## Reed-Solomon playground

```python
import algebrapy as alg

rs = alg.ReedSolomonCode(2, [1, 1, 0, 1], 3)   # GF(2^3), RS(7,3,5)
F = rs.field()
msg = [F.elem([1]), F.elem([0, 1]), F.elem([1, 1])]

codeword = rs.encode_systematic(msg)
print("codeword =", codeword)
print("syndromes =", rs.syndromes(codeword))

received = codeword.copy()
received[1] = received[1] + F.elem([1])
received[5] = received[5] + F.elem([1, 1])
print("decoded message =", rs.decode_message(received))
```

There is also a runnable example in [`play/rs.py`](/Users/froran/Projects/gautelis/algebra/play/rs.py).

## Setting up a Python environment

```terminaloutput
➜  python -m venv .venv
➜  source .venv/bin/activate
(.venv) ➜  pip install maturin
Collecting maturin
  Using cached maturin-1.12.4-py3-none-macosx_10_12_x86_64.macosx_11_0_arm64.macosx_10_12_universal2.whl.metadata (16 kB)
Using cached maturin-1.12.4-py3-none-macosx_10_12_x86_64.macosx_11_0_arm64.macosx_10_12_universal2.whl (18.9 MB)
Installing collected packages: maturin
Successfully installed maturin-1.12.4

```

## Building the Algebra backend
```terminaloutput
(.venv) ➜  (cd algebrapy && maturin develop)
🔗 Found pyo3 bindings
🐍 Found CPython 3.12 at /Users/froran/Projects/gautelis/algebra/.venv/bin/python
   Compiling pyo3-build-config v0.28.2
   Compiling pyo3-ffi v0.28.2
   Compiling pyo3-macros-backend v0.28.2
   Compiling pyo3 v0.28.2
   Compiling pyo3-macros v0.28.2
   Compiling algebrapy v0.1.0 (/Users/froran/Projects/gautelis/algebra/algebrapy)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 4.59s
📦 Built wheel for CPython 3.12 to /var/folders/ld/gmfpcx6n12j8w922_4_sbcvc0000gn/T/.tmpTO7MXa/algebrapy-0.1.0-cp312-cp312-macosx_11_0_arm64.whl
✏️ Setting installed package as editable
🛠 Installed algebrapy-0.1.0
```

For local Python smoke tests after `maturin develop`, use the virtualenv interpreter directly:

```terminaloutput
(.venv) ➜  .venv/bin/python -c "import algebrapy; print(algebrapy.Fp(7).elem(3))"
3 (mod 7)
```

On this repository, `python3` may still resolve to the system interpreter instead of `.venv/bin/python`, so `import algebrapy` can fail even though the editable install succeeded.

## Play with the Algebra package
```terminaloutput
(.venv) ➜  cd play
(.venv) ➜  python ghi.py
a = 3 (mod 7) b = 5 (mod 7)
a+b = 1 (mod 7)
a*b = 1 (mod 7)
inv(a) = 5 (mod 7)
a^6 = 1 (mod 7)
order(a) = 6
Fermat check OK

(.venv) ➜  cat rst.py 
import algebrapy as alg

K = alg.Fq(2, [1,1,0,1])      # 1 + x + x^3
a = K.elem([1,0,1])           # 1 + x^2
b = K.elem([0,1])             # x

print("a:", a)
print("b:", b)
print("a+b:", a + b)
print("a*b:", a * b)
print("b**7:", b**7)          # should be 1 in GF(2^3)
print("a/b:", a / b)
print("a**-1:", a**-1)
print("3 + a:", 3 + a)
print("2 * b:", 2 * b)        # note: in F2, 2 ≡ 0 => gives 0

(.venv) ➜  python rst.py 
a: 1 + x^2 in GF(2^3)
b: x in GF(2^3)
a+b: 1 + x + x^2 in GF(2^3)
a*b: 1 in GF(2^3)
b**7: 1 in GF(2^3)
a/b: 1 + x + x^2 in GF(2^3)
a**-1: x in GF(2^3)
3 + a: x^2 in GF(2^3)
2 * b: 0 in GF(2^3)

```

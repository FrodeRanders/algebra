# Algebra for Python

This is algebrapy - a Python package for abstract algebra implemented in Rust (using PyO3).

## What it provides

1. Finite Fields
- Fp / FpElem - Prime fields GF(p) where p is prime
- Fq / FqElem - Extension fields GF(p^k) constructed as polynomial quotient rings Fpx/(f) where f is irreducible
2. Finite Rings
- Zn / ZnElem - Integer residue rings Z/nZ with unit detection, zero-divisor checks, and inversion for units
3. Permutation Groups
- Perm - Permutations (bijections) with composition, inversion, exponentiation
- Sn - Symmetric groups S_n with element enumeration and subgroup generation
4. Arithmetic Utilities
- Extended GCD algorithm
- Prime number checking
- Modular inverse
5. Coding Theory
- BinaryBchCode - primitive narrow-sense binary BCH generator construction, encoding, and syndrome checks
- ReedSolomonCode - primitive narrow-sense Reed-Solomon codes over GF(p^m)

For a method-by-method catalog of the algebraic APIs, see [`docs/group-ring-field-catalog.md`](docs/group-ring-field-catalog.md).

## Example usage

Coefficient convention for `PolyFp`, `Fq`, and related coding constructors:
- coefficient lists are written low degree to high degree
- `[a0, a1, a2, a3]` means `a0 + a1*x + a2*x^2 + a3*x^3`

```python
import algebrapy as alg

# Prime field GF(7)
F = alg.Fp(7)
a = F.elem(3)
b = F.elem(5)

print("GF(7)")
print("elements =", F.elements())
print("a =", a, "b =", b)
print("a + b =", a + b)
print("a * b =", a * b)
print("a / b =", a / b)
print("a^6 =", a**6)          # Fermat: a^(p-1) = 1 for a != 0
print("order(a) =", a.mul_order())
print()

# Extension field GF(2^3) = F2[x] / (1 + x^2 + x^3)
# coefficients are low -> high, so [1, 0, 1, 1] means 1 + x^2 + x^3
K = alg.Fq(2, [1, 0, 1, 1])
x = K.elem([0, 1])            # x
c = K.elem([1, 0, 1])         # 1 + x^2

print("GF(2^3)")
print("x =", x)
print("c =", c)
print("x + c =", x + c)
print("x * c =", x * c)
print("x^7 =", x**7)
print("c^-1 =", c**-1)
print()

# Residue ring Z/12Z
R = alg.Zn(12)
u = R.elem(5)
z = R.elem(6)
d = R.elem(4)
e = R.elem(3)

print("Z/12Z")
print("modulus =", R.modulus())
print("elements =", R.elements())
print("u =", u, "z =", z)
print("u+z =", u + z)
print("u*z =", u * z)
print("u^3 =", u ** 3)
print("units =", R.units())
print("zero divisors =", R.zero_divisors())
print("R is integral domain?", R.is_integral_domain())
print("u is unit?", u.is_unit())
print("z is unit?", z.is_unit())
print("u is zero divisor?", u.is_zero_divisor())
print("z is zero divisor?", z.is_zero_divisor())
print("d =", d, "e =", e)
print("d*e =", d * e, "so both are zero divisors")
print("inv(u) =", u.inv())

try:
    print("inv(z) =", z.inv())
except Exception as exc:
    print("inv(z) failed:", exc)

print("3 + u =", 3 + u)
print("3 * u =", 3 * u)
print("3 / u =", 3 / u)
```

which yields the output:
```terminaloutput
GF(7)
elements = [0 (mod 7), 1 (mod 7), 2 (mod 7), 3 (mod 7), 4 (mod 7), 5 (mod 7), 6 (mod 7)]
a = 3 (mod 7) b = 5 (mod 7)
a + b = 1 (mod 7)
a * b = 1 (mod 7)
a / b = 2 (mod 7)
a^6 = 1 (mod 7)
order(a) = 6

GF(2^3)
x = x in GF(2^3)
c = 1 + x^2 in GF(2^3)
x + c = 1 + x + x^2 in GF(2^3)
x * c = 1 + x + x^2 in GF(2^3)
x^7 = 1 in GF(2^3)
c^-1 = 1 + x + x^2 in GF(2^3)

Z/12Z
modulus = 12
elements = [0 (mod 12), 1 (mod 12), 2 (mod 12), 3 (mod 12), 4 (mod 12), 5 (mod 12), 6 (mod 12), 7 (mod 12), 8 (mod 12), 9 (mod 12), 10 (mod 12), 11 (mod 12)]
u = 5 (mod 12) z = 6 (mod 12)
u+z = 11 (mod 12)
u*z = 6 (mod 12)
u^3 = 5 (mod 12)
units = [1 (mod 12), 5 (mod 12), 7 (mod 12), 11 (mod 12)]
zero divisors = [2 (mod 12), 3 (mod 12), 4 (mod 12), 6 (mod 12), 8 (mod 12), 9 (mod 12), 10 (mod 12)]
R is integral domain? False
u is unit? True
z is unit? False
u is zero divisor? False
z is zero divisor? True
d = 4 (mod 12) e = 3 (mod 12)
d*e = 0 (mod 12) so both are zero divisors
inv(u) = 5 (mod 12)
inv(z) failed: element is not a unit in Zn
3 + u = 8 (mod 12)
3 * u = 3 (mod 12)
3 / u = 3 (mod 12)
```

The implementation supports addition, subtraction, multiplication, and exponentiation across the algebraic structures it exposes. Division and negative exponents are available when an element is invertible, which for `Zn` means the element is a unit. For `Zn`, zero-divisor checks are exposed explicitly so you can experiment with rings that are not integral domains. Built with maturin for easy Python integration.

There is also a runnable example in [`play/zn_ring_basics.py`](play/zn_ring_basics.py).

For the relationship between residue rings and permutation groups, see [`play/zn_permutation_actions.py`](play/zn_permutation_actions.py). It compares prime and composite moduli using translations `x -> x + b`, unit actions `x -> a*x`, affine maps `x -> a*x + b`, cycle decompositions, and orbit structure.

## BCH playground

```python
import algebrapy as alg

# Primitive binary BCH / Hamming code of length 7 over GF(2^3)
# polynomial coefficients are low -> high
code = alg.BinaryBchCode(3, [1, 1, 0, 1], 3)   # 1 + x + x^3
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

There is also a fuller runnable example in [`play/bch_code_demo.py`](play/bch_code_demo.py) that exercises:
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

# modulus coefficients are low -> high, so [1, 1, 0, 1] means 1 + x + x^3
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

There is also a runnable example in [`play/reed_solomon_demo.py`](play/reed_solomon_demo.py).

## CD-style scratch demo

Compact discs use a more elaborate Reed-Solomon-based pipeline called CIRC, but the essential idea can be demonstrated with the current codebase: interleaving spreads a localized burst error across several codewords so each individual Reed-Solomon decoder only sees a small number of symbol errors.

The toy demo in [`play/reed_solomon_cd_scratch_demo.py`](play/reed_solomon_cd_scratch_demo.py) shows exactly that:
- without interleaving, a short burst damages too many symbols in one codeword and decoding fails
- with interleaving, the same burst is distributed across several codewords and each one decodes successfully

This is intentionally not a literal CD implementation. Real CD audio/data protection uses:
- two Reed-Solomon layers rather than one
- byte-oriented symbols over `GF(2^8)`
- shortened codes and deeper interleaving
- additional framing and concealment machinery around the codes

So the demo is best read as "why Reed-Solomon plus interleaving helps with scratches", not as a full simulation of the exact CD standard.

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
🐍 Found CPython 3.12 at .venv/bin/python
   Compiling pyo3-build-config v0.28.2
   Compiling pyo3-ffi v0.28.2
   Compiling pyo3-macros-backend v0.28.2
   Compiling pyo3 v0.28.2
   Compiling pyo3-macros v0.28.2
   Compiling algebrapy v0.1.0 (./algebrapy)
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

Example scripts in `play/`:
- `fp_field_methods.py` - prime-field operations through the explicit `Fp` API
- `fp_operator_overloads.py` - prime-field arithmetic via Python operators
- `fp_pow_protocol.py` - exponentiation and `pow(...)` on field elements
- `fp_multiplicative_orders.py` - multiplicative orders and primitive elements in `Fp*`
- `fp_quadratic_residue_table.py` - quadratic residue checks and square roots over `Fp`
- `fp_quadratic_residues_and_dlog.py` - quadratic residues together with discrete logarithms
- `fp_cayley_table.py` - multiplication table for a small prime field
- `fq_extension_construction.py` - building `GF(p^k)` from an irreducible polynomial
- `fq_operator_overloads.py` - extension-field arithmetic via Python operators
- `fq_power_cycle.py` - powers of a distinguished element in an extension field
- `fq_trace_and_norm.py` - trace and norm in `GF(2^3)`
- `fq_trace_norm_over_f9.py` - trace, norm, and primitive elements in `GF(3^2)`
- `fq_primitive_elements.py` - primitive elements in a small extension field
- `zn_ring_basics.py` - arithmetic, units, zero divisors, and inversion in `Z/nZ`
- `zn_permutation_actions.py` - translations, unit actions, affine permutations, and orbit structure from `Z/nZ` into `S_n`
- `sn_subgroup_generation.py` - subgroup generation in a symmetric group
- `bch_code_demo.py` - BCH encoding, decoding, shortening, and code parameters
- `reed_solomon_demo.py` - Reed-Solomon encoding, syndromes, and decoding

```terminaloutput
(.venv) ➜  cd play
(.venv) ➜  python fp_field_methods.py
field modulus = 7
elements = [0 (mod 7), 1 (mod 7), 2 (mod 7), 3 (mod 7), 4 (mod 7), 5 (mod 7), 6 (mod 7)]
nonzero elements = [1 (mod 7), 2 (mod 7), 3 (mod 7), 4 (mod 7), 5 (mod 7), 6 (mod 7)]
a = 3 (mod 7) b = 5 (mod 7)
a+b = 1 (mod 7)
a*b = 1 (mod 7)
inv(a) = 5 (mod 7)
a^6 = 1 (mod 7)
order(a) = 6
Fermat check OK

(.venv) ➜  cat fq_operator_overloads.py
import algebrapy as alg

K = alg.Fq(2, [1, 1, 0, 1])  # coefficients are low -> high: 1 + x + x^3
a = K.elem([1, 0, 1])  # 1 + x^2
b = K.elem([0, 1])  # x

print("a:", a)
print("b:", b)
print("a+b:", a + b)
print("a-b:", a - b)
print("a*b:", a * b)
print("b**7:", b**7)
print("a/b:", a / b)
print("a**-1:", a**-1)
print("3 + a:", 3 + a)
print("2 * b:", 2 * b)

(.venv) ➜  python fq_operator_overloads.py
a: 1 + x^2 in GF(2^3)
b: x in GF(2^3)
a+b: 1 + x + x^2 in GF(2^3)
a-b: 1 + x + x^2 in GF(2^3)
a*b: 1 in GF(2^3)
b**7: 1 in GF(2^3)
a/b: 1 + x + x^2 in GF(2^3)
a**-1: x in GF(2^3)
3 + a: x^2 in GF(2^3)
2 * b: 0 in GF(2^3)

```

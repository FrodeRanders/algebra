# Scaling Plan for Larger Algebraic Objects

This document is a roadmap for growing `algebrapy` beyond its current
small-explicit-object model. It is intentionally a planning document, not a
commitment to implement every item.

The current project is strongest as a concrete finite algebra playground:
objects are enumerable, operations are inspectable, and group actions are often
materialized as actual permutations. That is excellent for learning and for
small examples, but larger algebraic structures need a different default:
compact structural representations, lazy algorithms, and explicit bounds around
enumeration.

## Current Scaling Limits

The main limitations are not Rust or Python overhead. They are mathematical
representation choices:

- Many APIs assume a finite object can be fully enumerated.
- Arithmetic actions such as `x -> a*x` are often converted into full
  permutation image lists.
- `PermSubgroup` stores every subgroup element explicitly.
- Some algorithms still use simple trial division or repeated multiplication.
- Polynomial arithmetic is dense and straightforward, with no specialized
  factorization or modular-composition machinery.

This is appropriate for examples such as `GF(2^3)`, `Z/12Z`, `S4`, BCH(15,7,5),
and RS(7,3,5). It does not scale to objects such as `GF(2^128)`, large
Reed-Solomon fields, or permutation groups with huge generated order.

## Design Principle

The central architectural shift should be:

> A finite algebraic object has compact structure; enumeration is an optional
> operation with explicit bounds.

For example, `GF(2^128)` should support arithmetic, inverses, trace, norm,
minimal-polynomial style operations, and coding-theory use cases without ever
calling `elements()`. Asking for `elements()` should remain possible only under
explicit size limits.

## API Shape

### Keep Small Explicit APIs

The current explicit APIs are valuable and should remain:

- `elements()`
- `add_perm()`, `mul_perm()`, `affine_perm()`
- subgroup enumeration
- cycle notation and orbit decomposition by materialized permutations
- codeword enumeration and weight distributions for small codes

These should be documented as small-object APIs and guarded consistently by
`max_size` or similarly named bounds.

### Add Non-Enumerating APIs

Large-object support should add methods that compute from structure:

- `order()` without enumeration where mathematically known.
- `is_unit()` by gcd or field nonzeroness.
- `mul_order()` using factorization of group order.
- `is_primitive()` using prime divisors of `q - 1`.
- `trace()` and `norm()` without enumerating field elements.
- `minimal_polynomial()` without enumerating all polynomials.
- action objects that can apply to a point without materializing a full
  permutation.

### Make Enumeration Explicit

Potentially large operations should advertise their cost in the name or
signature:

- `elements(max_size=...)`
- `as_permutation(max_size=...)`
- `orbits_explicit(max_size=...)`
- `codewords(max_size=...)`
- `weight_distribution(max_codewords=...)`

For APIs that currently use defaults, prefer conservative defaults and helpful
errors that suggest the non-enumerating alternative.

## Finite Fields

### Irreducibility

`PolyFp.is_irreducible()` now uses Rabin's irreducibility test. That is a good
first step. Further improvements would be:

- Use modular composition for faster Frobenius computations.
- Avoid repeated polynomial exponentiation when checking many candidate
  polynomials of the same degree.
- Add separate tests for square-free, irreducible, and primitive polynomials.

### Primitive Elements and Multiplicative Order

The current `Fq.mul_order()` is simple repeated multiplication up to `q - 1`.
For large fields, use factorization of `q - 1`:

1. Factor `q - 1`.
2. For an element `a`, start with `order = q - 1`.
3. For each prime divisor `r`, repeatedly test whether
   `a^(order / r) == 1`.
4. Divide `order` by `r` when the test succeeds.

Primitive-element search should similarly test whether
`a^((q - 1) / r) != 1` for every prime divisor `r` of `q - 1`.

This requires integer factorization utilities. Trial division may be fine at
first, but Pollard rho or a dependency should be considered for larger orders.

### Representation

The current polynomial-basis representation is reasonable and should remain:

- coefficients in `F_p`
- reduction modulo a monic irreducible polynomial
- zero as an empty coefficient vector

Possible upgrades:

- fixed-width representations for binary fields, especially `GF(2^m)`;
- specialized carryless multiplication for `p = 2`;
- cached modulus metadata such as degree, reciprocal, and factorization of
  `q - 1`;
- faster equality and hashing for high-degree field elements.

## Polynomial Arithmetic

The current dense-vector implementation is clear and adequate for small
degrees. Larger work would benefit from:

- faster gcd algorithms;
- square-free factorization;
- distinct-degree factorization;
- equal-degree factorization, such as Cantor-Zassenhaus;
- modular composition;
- sparse representation for high-degree low-density polynomials;
- optional specialized binary-polynomial backend.

This area is a good candidate for an internal trait-like split:

- simple dense implementation for clarity and small sizes;
- optimized implementation selected internally for larger degrees or `p = 2`.

## Residue Rings

`Zn` already has a compact representation, so it scales better than explicit
fields in some ways. The main upgrades are number-theoretic:

- factor `n`;
- compute Euler phi and Carmichael lambda;
- expose Chinese remainder decomposition when `n` is factored;
- compute unit group structure for supported moduli;
- avoid enumerating units for large `n`;
- provide lazy unit actions rather than full permutation representations.

For example, `Z/nZ` can answer whether an element is a unit using `gcd(a, n)`
without enumerating anything. But `unit_group_perms()` should remain
small-object-only.

## Arithmetic Actions

The current finite-action viewpoint is central to the project and should not be
removed. It should be generalized.

### Current Model

Arithmetic action -> materialized `Perm`.

This works for `GF(7)`, `GF(2^3)`, and `Z/12Z`.

### Larger Model

Arithmetic action -> lazy action object.

Possible type:

```text
FiniteAction {
    domain metadata,
    apply(point),
    compose(other),
    inverse() when available,
    maybe_as_perm(max_size),
}
```

This allows APIs such as:

- `field.mul_action(a)`
- `ring.affine_action(a, b)`
- `action.apply(x)`
- `action.as_perm(max_size=4096)`

For small domains, the action can still be converted to `Perm`. For large
domains, it remains a callable algebraic map.

## Permutation Groups

The current `PermSubgroup` is explicit: it stores every element. That is the
main barrier to larger permutation groups.

A scalable permutation-group layer would need Schreier-Sims machinery:

- base and strong generating set;
- stabilizer chains;
- membership testing without full enumeration;
- group order computation;
- random element generation;
- orbit computations from generators;
- point stabilizers without enumerating all group elements.

This is a significant project. There are two practical paths:

1. Implement a modest Schreier-Sims layer in Rust for educational and medium
   scale use.
2. Bind to or interoperate with a mature computational algebra system for large
   permutation groups.

The current explicit subgroup type should remain as `ExplicitPermSubgroup` or
continue as `PermSubgroup` with clear size limits. A future non-enumerating
group type could be introduced alongside it.

## Coding Theory

BCH and Reed-Solomon code support already depends on finite fields and
polynomial arithmetic, so it benefits directly from field and polynomial
upgrades.

Important scaling steps:

- avoid codeword enumeration except under explicit bounds;
- keep syndrome, encoding, and decoding non-enumerating;
- use faster finite-field arithmetic for common fields such as `GF(2^8)`;
- separate educational matrix construction from production encoding paths;
- add benchmarks for common code sizes;
- avoid cloning full field objects unnecessarily in inner loops.

For Reed-Solomon codes, larger support particularly wants efficient symbol
arithmetic over `GF(2^8)` and possibly `GF(2^16)`.

## Error Handling and Documentation

Large-object support should make cost visible:

- methods that enumerate should say so in documentation;
- errors should explain which size bound was exceeded;
- non-enumerating alternatives should be suggested in error messages;
- examples should distinguish "small explicit demo" from "large structural use".

Example error:

```text
field too large to enumerate as a permutation; use mul_action(a) or increase max_size
```

## Performance Work

Before optimizing aggressively, add benchmarks for:

- polynomial multiplication and division;
- polynomial gcd;
- Rabin irreducibility testing;
- `Fq` multiplication, inversion, trace, and norm;
- primitive-element search;
- BCH generator construction;
- Reed-Solomon encode/decode;
- permutation subgroup generation on known small groups.

Use those benchmarks to guide implementation rather than guessing.

## Suggested Implementation Phases

### Phase 1: Make Existing Boundaries Explicit

- Audit methods that enumerate.
- Add or standardize `max_size` parameters.
- Improve error messages for oversized enumeration.
- Document small-object APIs in the method catalog.
- Add tests for refusing oversized operations.

### Phase 2: Finish Field and Polynomial Algorithm Upgrades

- Keep Rabin's irreducibility test.
- Add integer factorization helpers.
- Replace `Fq.mul_order()` with factorization-based order computation.
- Replace primitive-element search with factorization-based tests.
- Add polynomial factorization primitives.

### Phase 3: Introduce Lazy Arithmetic Actions

- Add action objects for `Fp`, `Fq`, and `Zn`.
- Keep `add_perm`, `mul_perm`, and `affine_perm` as small-domain conversions.
- Add `as_perm(max_size)` to action objects.
- Add examples showing the same action used lazily and explicitly.

### Phase 4: Larger Coding-Theory Use

- Optimize common binary extension fields.
- Avoid unnecessary field enumeration in coding APIs.
- Add benchmarks and examples for larger RS/BCH parameters.
- Separate matrix and codeword enumeration helpers from normal encode/decode
  workflows.

### Phase 5: Non-Enumerating Permutation Groups

- Decide whether to implement Schreier-Sims or integrate with an existing
  system.
- Add a non-enumerating generated permutation group type.
- Keep explicit subgroup enumeration for small examples.
- Rework Sylow and normality APIs around clear capability limits.

## Open Decisions

- Should `PolyFp` continue accepting composite `p` for general polynomial rings,
  while field-specific methods require prime `p`?
- Should optimized binary-field arithmetic live inside `Fq`, or should there be
  a separate `GF2m` type?
- Should large permutation-group support be native Rust, optional dependency, or
  delegated to an external system?
- How much of the API should prioritize pedagogy versus performance?
- Should action objects become a shared abstraction used by fields, rings, and
  permutation groups?

## Near-Term Candidates

The most practical next changes are:

1. Add integer factorization helpers.
2. Replace `Fq.mul_order()` with factorization-based order computation.
3. Replace `Fq.primitive_elements()` internals with primitive-element tests.
4. Add `max_size` tests around all enumerating APIs.
5. Prototype lazy arithmetic actions for `Zn` and `Fp`.

Those steps preserve the current educational feel while moving the codebase
toward larger finite algebra.

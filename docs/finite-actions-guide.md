# Finite Actions Guide

This guide collects the action-oriented parts of `algebrapy` into one place.

The basic idea is simple:
- algebraic operations such as `x -> x + b` and `x -> a*x` can be viewed as permutations of the underlying finite set
- once you have those permutations, the group API can tell you about cycles, orbits, subgroup sizes, transitivity, and stabilizers

That gives a concrete way to compare fields, extension fields, and rings with zero divisors.

## Core viewpoint

For a finite set `X`, a permutation action is just a way to let algebraic data move points of `X` around.

In this codebase:
- `Perm` represents a single permutation
- `Sn` provides subgroup and action utilities on `{0, ..., n-1}`
- `Fp`, `Fq`, and `Zn` can build permutations from arithmetic maps

At this point, the arithmetic side is fairly parallel:
- `Fp` exposes multiplicative-action summaries directly
- small `Fq` exposes the same summaries, with enumeration bounds
- `Zn` exposes the corresponding unit-action summaries

The most important arithmetic maps are:
- translation: `x -> x + b`
- multiplication: `x -> a*x`
- affine maps: `x -> a*x + b`

## Why fields and rings behave differently

### Prime fields `GF(p)`

In a field, every nonzero element is invertible, so every nonzero `a` defines a permutation

`x -> a*x`

on the full underlying set.

Consequences:
- multiplication by a primitive element cycles through all nonzero elements
- the nonzero part of the field often appears as one large orbit under multiplicative action
- affine maps behave especially cleanly

Relevant API:
- `Fp.add_perm`
- `Fp.mul_perm`
- `Fp.affine_perm`
- `Fp.mul_order`

Runnable example:
- [`play/fp_permutation_actions.py`](play/fp_permutation_actions.py)

### Extension fields `GF(p^k)`

The same picture holds for small extension fields, with one practical difference: to turn arithmetic into an explicit permutation, the implementation must enumerate the whole field.

That is why the `Fq` action methods have bounded variants:
- `add_perm_with_limit`
- `mul_perm_with_limit`
- `affine_perm_with_limit`

For small fields such as `GF(2^3)`, this works well and makes primitive-element actions very concrete.

Relevant API:
- `Fq.add_perm`
- `Fq.mul_perm`
- `Fq.affine_perm`
- `Fq.primitive_elements`
- `Fq.mul_order`

Runnable example:
- [`play/fq_permutation_actions.py`](play/fq_permutation_actions.py)

### Residue rings `Z/nZ`

This is where the action viewpoint becomes especially educational.

In `Z/nZ`:
- `x -> x + b` is always a permutation
- `x -> a*x` is a permutation if and only if `a` is a unit
- zero divisors fail to give permutations

So permutation behavior becomes a direct test for the difference between units and zero divisors.

Relevant API:
- `Zn.units`
- `Zn.zero_divisors`
- `Zn.is_integral_domain`
- `Zn.add_perm`
- `Zn.mul_perm`
- `Zn.affine_perm`
- `Zn.unit_group_perms`
- `Zn.unit_action_subgroup_size`
- `Zn.unit_action_orbits`
- `Zn.unit_action_stabilizer`
- `Zn.unit_action_stabilizer_size`
- `Zn.is_unit_action_transitive`

Runnable examples:
- [`play/zn_ring_basics.py`](play/zn_ring_basics.py)
- [`play/zn_permutation_actions.py`](play/zn_permutation_actions.py)

## Reading the permutation data

Once you have a permutation, the most useful inspection methods are:
- `Perm.cycle_notation()`
- `Perm.cycle_type()`
- `Perm.order()`

Examples:
- a translation in `Z/7Z` often appears as a single 7-cycle
- multiplication by a primitive element in `GF(7)` fixes `0` and permutes the six nonzero elements in one 6-cycle
- multiplication by a unit in `Z/12Z` usually breaks into several smaller cycles

That already shows the algebra.

## Reading the group action data

Given one or more generators, `Sn` can tell you how the generated subgroup acts:

- `orbit(point, gens)`:
  where one point can move
- `orbits(gens)`:
  how the whole set splits into invariant pieces
- `subgroup_size(gens)`:
  how large the generated action group is
- `is_transitive(gens)`:
  whether there is only one orbit
- `stabilizer(point, gens)` and `stabilizer_size(point, gens)`:
  which elements of the subgroup fix a chosen point

This is the standard orbit-stabilizer toolkit in a concrete computational form.

## Beyond arithmetic: conjugation in `S_n`

The same permutation tools are also useful without any ring or field in the picture.

Inside `S_n`, permutations act on each other by conjugation:

`g * p * g^-1`

This does not change cycle type, which is why cycle structure organizes conjugacy classes in symmetric groups.

Relevant API:
- `Perm.conjugate_by`
- `Sn.conjugacy_class`
- `Sn.conjugacy_class_size`

Runnable example:
- [`play/sn_conjugacy_demo.py`](play/sn_conjugacy_demo.py)

## A good comparison to keep in mind

### `GF(7)`

Take `g = 3`.

Then:
- `g` has multiplicative order `6`
- `x -> 3x` fixes `0`
- on nonzero elements, the orbit of `1` is `[1, 2, 3, 4, 5, 6]`

So the multiplicative action is cyclic on the nonzero part.

### `Z/12Z`

Take the full unit group `{1, 5, 7, 11}`.

Then the action on residues splits into several orbits:
- `[0]`
- `[1, 5, 7, 11]`
- `[2, 10]`
- `[3, 9]`
- `[4, 8]`
- `[6]`

The stabilizers also differ:
- `1` has stabilizer size `1`
- `6` has stabilizer size `4`

That fragmentation reflects the presence of zero divisors and non-units.

## Suggested study path

If you want to learn the codebase in a mathematically sensible order:

1. Start with [`play/zn_ring_basics.py`](play/zn_ring_basics.py)
   Focus on units, zero divisors, and why inverses fail.
2. Move to [`play/fp_permutation_actions.py`](play/fp_permutation_actions.py)
   Focus on one big nonzero orbit and primitive elements.
3. Compare with [`play/zn_permutation_actions.py`](play/zn_permutation_actions.py)
   Focus on orbit splitting in `Z/12Z`.
4. Then read [`play/fq_permutation_actions.py`](play/fq_permutation_actions.py)
   Focus on how the same picture extends to `GF(2^3)`.
5. Finally read [`play/sn_conjugacy_demo.py`](play/sn_conjugacy_demo.py)
   Focus on what the same permutation toolkit says about groups acting on themselves by conjugation.

## Where to look next

For API lookup:
- [`docs/group-ring-field-catalog.md`](docs/group-ring-field-catalog.md)

For runnable demonstrations:
- [`play/fp_permutation_actions.py`](play/fp_permutation_actions.py)
- [`play/fq_permutation_actions.py`](play/fq_permutation_actions.py)
- [`play/zn_permutation_actions.py`](play/zn_permutation_actions.py)
- [`play/sn_conjugacy_demo.py`](play/sn_conjugacy_demo.py)

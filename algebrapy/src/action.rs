use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::arith::egcd::inv_mod_i128;
use crate::field::fq::{Fq, FqElem};
use crate::group::perm::Perm;

/// The concrete arithmetic map behind a `FiniteAction`.
///
/// Every variant is stored in normal form as an affine map `x -> a*x + b`,
/// where `a` is invertible in the domain.
#[derive(Clone, Debug)]
enum ActionKind {
    Fp { p: u64, a: u64, b: u64 },
    Fq { field: Fq, a: FqElem, b: FqElem },
    Zn { n: u64, a: u64, b: u64 },
}

fn classify(a_is_one: bool, b_is_zero: bool) -> &'static str {
    match (a_is_one, b_is_zero) {
        (true, true) => "identity",
        (true, false) => "translation",
        (false, true) => "multiplication",
        (false, false) => "affine",
    }
}

/// A lazy permutation action `x -> a*x + b` on a finite field or residue ring.
///
/// Points are encoded as `0, 1, ..., size-1`. The action applies arithmetic
/// directly to a point, so it never has to enumerate or store the domain.
/// Use `as_perm(max_size=...)` to materialize it as a `Perm` when the domain
/// is small enough.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct FiniteAction {
    kind: ActionKind,
}

impl FiniteAction {
    pub(crate) fn fp(p: u64, a: u64, b: u64) -> Self {
        Self {
            kind: ActionKind::Fp {
                p,
                a: a % p,
                b: b % p,
            },
        }
    }

    pub(crate) fn fq(field: Fq, a: FqElem, b: FqElem) -> Self {
        Self {
            kind: ActionKind::Fq { field, a, b },
        }
    }

    pub(crate) fn zn(n: u64, a: u64, b: u64) -> Self {
        Self {
            kind: ActionKind::Zn {
                n,
                a: a % n,
                b: b % n,
            },
        }
    }

    fn domain_size(&self) -> u64 {
        match &self.kind {
            ActionKind::Fp { p, .. } => *p,
            ActionKind::Fq { field, .. } => field.size(),
            ActionKind::Zn { n, .. } => *n,
        }
    }

    fn kind_name(&self) -> &'static str {
        match &self.kind {
            ActionKind::Fp { a, b, .. } | ActionKind::Zn { a, b, .. } => classify(*a == 1, *b == 0),
            ActionKind::Fq { a, b, .. } => classify(a.coeffs() == vec![1], b.is_zero()),
        }
    }

    fn apply_internal(&self, point: u64) -> PyResult<u64> {
        match &self.kind {
            ActionKind::Fp { p, a, b } => {
                if point >= *p {
                    return Err(PyValueError::new_err("point out of range"));
                }
                Ok((((*a as u128) * (point as u128) + (*b as u128)) % (*p as u128)) as u64)
            }
            ActionKind::Zn { n, a, b } => {
                if point >= *n {
                    return Err(PyValueError::new_err("point out of range"));
                }
                Ok((((*a as u128) * (point as u128) + (*b as u128)) % (*n as u128)) as u64)
            }
            ActionKind::Fq { field, a, b } => {
                if point >= field.size() {
                    return Err(PyValueError::new_err("point out of range"));
                }
                let x = field.point_to_elem(point);
                let ax = field.mul(a, &x)?;
                let y = field.add(&ax, b)?;
                Ok(field.elem_to_point(&y))
            }
        }
    }
}

#[pymethods]
impl FiniteAction {
    /// Return the number of points in the domain.
    pub fn size(&self) -> u64 {
        self.domain_size()
    }

    /// Return the action type: `identity`, `translation`, `multiplication`, or `affine`.
    pub fn kind(&self) -> &'static str {
        self.kind_name()
    }

    /// Apply the action to `point` without enumerating the domain.
    pub fn apply(&self, point: u64) -> PyResult<u64> {
        self.apply_internal(point)
    }

    /// Return the cycle of `point` under repeated application.
    pub fn cycle(&self, point: u64) -> PyResult<Vec<u64>> {
        let size = self.domain_size();
        if point >= size {
            return Err(PyValueError::new_err("point out of range"));
        }
        let mut out = vec![point];
        let mut current = self.apply_internal(point)?;
        while current != point {
            if out.len() as u64 > size {
                return Err(PyValueError::new_err("action is not a permutation"));
            }
            out.push(current);
            current = self.apply_internal(current)?;
        }
        Ok(out)
    }

    /// Return the composition `self ∘ other` (apply `other`, then `self`).
    pub fn compose(&self, other: &FiniteAction) -> PyResult<FiniteAction> {
        match (&self.kind, &other.kind) {
            (
                ActionKind::Fp { p, a: a1, b: b1 },
                ActionKind::Fp {
                    p: p2,
                    a: a2,
                    b: b2,
                },
            ) => {
                if p != p2 {
                    return Err(PyValueError::new_err("actions have different domains"));
                }
                let p = *p;
                let a = ((*a1 as u128) * (*a2 as u128) % (p as u128)) as u64;
                let b = (((*a1 as u128) * (*b2 as u128) + (*b1 as u128)) % (p as u128)) as u64;
                Ok(Self::fp(p, a, b))
            }
            (
                ActionKind::Zn { n, a: a1, b: b1 },
                ActionKind::Zn {
                    n: n2,
                    a: a2,
                    b: b2,
                },
            ) => {
                if n != n2 {
                    return Err(PyValueError::new_err("actions have different domains"));
                }
                let n = *n;
                let a = ((*a1 as u128) * (*a2 as u128) % (n as u128)) as u64;
                let b = (((*a1 as u128) * (*b2 as u128) + (*b1 as u128)) % (n as u128)) as u64;
                Ok(Self::zn(n, a, b))
            }
            (
                ActionKind::Fq {
                    field,
                    a: a1,
                    b: b1,
                },
                ActionKind::Fq {
                    field: field2,
                    a: a2,
                    b: b2,
                },
            ) => {
                if !field.same_field(field2) {
                    return Err(PyValueError::new_err("actions have different domains"));
                }
                let a = field.mul(a1, a2)?;
                let b = field.add(&field.mul(a1, b2)?, b1)?;
                Ok(Self::fq(field.clone(), a, b))
            }
            _ => Err(PyValueError::new_err("actions have different domains")),
        }
    }

    /// Return the inverse action.
    pub fn inverse(&self) -> PyResult<FiniteAction> {
        match &self.kind {
            ActionKind::Fp { p, a, b } => {
                let inv = inv_mod_i128(*a as i128, *p as i128)
                    .ok_or_else(|| PyValueError::new_err("action is not invertible"))?
                    as u64;
                let p = *p;
                let b_new = ((p as u128 - ((inv as u128) * (*b as u128) % (p as u128)))
                    % (p as u128)) as u64;
                Ok(Self::fp(p, inv, b_new))
            }
            ActionKind::Zn { n, a, b } => {
                let inv = inv_mod_i128(*a as i128, *n as i128)
                    .ok_or_else(|| PyValueError::new_err("action is not invertible"))?
                    as u64;
                let n = *n;
                let b_new = ((n as u128 - ((inv as u128) * (*b as u128) % (n as u128)))
                    % (n as u128)) as u64;
                Ok(Self::zn(n, inv, b_new))
            }
            ActionKind::Fq { field, a, b } => {
                let inv = field.inv(a)?;
                let b_new = field.neg(&field.mul(&inv, b)?)?;
                Ok(Self::fq(field.clone(), inv, b_new))
            }
        }
    }

    /// Materialize the action as a `Perm`, subject to `max_size`.
    #[pyo3(signature = (max_size=None))]
    pub fn as_perm(&self, max_size: Option<u64>) -> PyResult<Perm> {
        let max = max_size.unwrap_or(4096);
        let size = self.domain_size();
        if size > max {
            return Err(PyValueError::new_err(format!(
                "action domain of size {} exceeds max_size {}; use apply(point) or cycle(point), or increase max_size",
                size, max
            )));
        }
        let mut images: Vec<usize> = Vec::with_capacity(size as usize);
        for point in 0..size {
            images.push(self.apply_internal(point)? as usize);
        }
        Perm::new(size as usize, images)
    }

    pub fn __repr__(&self) -> String {
        match &self.kind {
            ActionKind::Fp { p, a, b } => {
                format!("FiniteAction(x -> {}*x + {} on GF({}))", a, b, p)
            }
            ActionKind::Zn { n, a, b } => {
                format!("FiniteAction(x -> {}*x + {} on Z/{}Z)", a, b, n)
            }
            ActionKind::Fq { field, a, b } => format!(
                "FiniteAction(x -> {:?}*x + {:?} on GF({}^{}))",
                a.coeffs(),
                b.coeffs(),
                field.p(),
                field.degree()
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::fp::Fp;
    use crate::ring::zn::Zn;

    #[test]
    fn fp_actions_apply_compose_and_invert() {
        let f = Fp::new(7).unwrap();
        let add = f.add_action(&f.elem(2)).unwrap();
        let mul = f.mul_action(&f.elem(3)).unwrap();
        let affine = f.affine_action(&f.elem(3), &f.elem(2)).unwrap();

        assert_eq!(add.kind(), "translation");
        assert_eq!(mul.kind(), "multiplication");
        assert_eq!(affine.kind(), "affine");
        assert_eq!(add.apply(6).unwrap(), 1);
        assert_eq!(mul.apply(5).unwrap(), 1);
        assert_eq!(affine.apply(2).unwrap(), 1);

        let composed = add.compose(&mul).unwrap();
        for x in 0..7 {
            assert_eq!(composed.apply(x).unwrap(), affine.apply(x).unwrap());
        }

        let inverse = affine.inverse().unwrap();
        for x in 0..7 {
            let y = affine.apply(x).unwrap();
            assert_eq!(inverse.apply(y).unwrap(), x);
        }
    }

    #[test]
    fn fp_action_matches_materialized_permutation() {
        let f = Fp::new(7).unwrap();
        let a = f.elem(3);
        let b = f.elem(2);

        assert_eq!(
            f.affine_action(&a, &b).unwrap().as_perm(Some(7)).unwrap(),
            f.affine_perm(&a, &b).unwrap()
        );
        assert_eq!(
            f.mul_action(&a).unwrap().as_perm(Some(7)).unwrap(),
            f.mul_perm(&a).unwrap()
        );
        assert_eq!(
            f.add_action(&b).unwrap().as_perm(Some(7)).unwrap(),
            f.add_perm(&b).unwrap()
        );
        assert!(f.add_action(&b).unwrap().as_perm(Some(3)).is_err());
    }

    #[test]
    fn zn_actions_require_units() {
        let z = Zn::new(12).unwrap();
        assert_eq!(z.mul_action(&z.elem(5)).unwrap().apply(7).unwrap(), 11);
        assert!(z.mul_action(&z.elem(6)).is_err());
        assert!(z.affine_action(&z.elem(6), &z.elem(1)).is_err());

        let affine = z.affine_action(&z.elem(5), &z.elem(1)).unwrap();
        assert_eq!(affine.apply(2).unwrap(), 11);
        let inverse = affine.inverse().unwrap();
        for x in 0..12 {
            let y = affine.apply(x).unwrap();
            assert_eq!(inverse.apply(y).unwrap(), x);
        }
    }

    #[test]
    fn fq_action_matches_materialized_permutation() {
        let field = Fq::new(2, vec![1, 1, 0, 0, 1]).unwrap();
        let a = field.elem(vec![1, 1]).unwrap();
        let b = field.elem(vec![1, 0, 1]).unwrap();

        assert_eq!(
            field
                .affine_action(&a, &b)
                .unwrap()
                .as_perm(Some(16))
                .unwrap(),
            field.affine_perm(&a, &b).unwrap()
        );
        assert_eq!(
            field.mul_action(&a).unwrap().as_perm(Some(16)).unwrap(),
            field.mul_perm(&a).unwrap()
        );

        let composed = field
            .add_action(&b)
            .unwrap()
            .compose(&field.mul_action(&a).unwrap())
            .unwrap();
        for x in 0..16 {
            assert_eq!(
                composed.apply(x).unwrap(),
                field.affine_action(&a, &b).unwrap().apply(x).unwrap()
            );
        }
    }

    #[test]
    fn large_fq_action_without_enumeration() {
        // GF(2^16) with x^16 + x^12 + x^3 + x + 1.
        let modulus: Vec<i128> = vec![1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1];
        let field = Fq::new(2, modulus).unwrap();
        assert_eq!(field.size(), 65536);

        let shift = field.add_action(&field.one()).unwrap();
        assert_eq!(shift.apply(0).unwrap(), 1);
        assert_eq!(shift.apply(1).unwrap(), 0);
        assert_eq!(shift.cycle(0).unwrap(), vec![0, 1]);

        // x is primitive for this modulus, so multiplication by x cycles
        // through all 2^16 - 1 nonzero elements.
        let x = field.elem(vec![0, 1]).unwrap();
        let mul = field.mul_action(&x).unwrap();
        let cycle = mul.cycle(1).unwrap();
        assert_eq!(cycle.len(), 65535);
        assert!(!cycle.contains(&0));
        assert!(mul.as_perm(Some(1024)).is_err());
    }
}

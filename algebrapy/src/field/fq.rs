use pyo3::exceptions::{PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::{Py, PyAny};

use super::poly_fp::PolyFp;
use crate::group::perm::Perm;

/// The extension field `GF(p^k)` represented as `F_p[x] / (f)`.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct Fq {
    p: u64,
    modulus: PolyFp, // irreducible polynomial
    k: usize,        // degree(modulus)
}

/// An element of a fixed extension field `GF(p^k)`.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct FqElem {
    p: u64,
    modulus_coeffs: Vec<u64>,
    coeffs: Vec<u64>, // reduced representative degree < k
}

fn trim(mut v: Vec<u64>) -> Vec<u64> {
    while let Some(&last) = v.last() {
        if last == 0 {
            v.pop();
        } else {
            break;
        }
    }
    v
}

fn poly_add(p: u64, a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut r = vec![0u64; n];
    for i in 0..n {
        let av = if i < a.len() { a[i] } else { 0 };
        let bv = if i < b.len() { b[i] } else { 0 };
        r[i] = (av + bv) % p;
    }
    trim(r)
}

fn poly_neg(p: u64, a: &[u64]) -> Vec<u64> {
    trim(a.iter().map(|&c| if c == 0 { 0 } else { p - c }).collect())
}

fn poly_mul(p: u64, a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            let add = ((a[i] as u128 * b[j] as u128) % (p as u128)) as u64;
            r[i + j] = (r[i + j] + add) % p;
        }
    }
    trim(r)
}

// Reduce poly a modulo modulus (given as coeffs, monic assumed is nice but not required)
fn poly_mod(p: u64, mut a: Vec<u64>, modulus: &[u64]) -> PyResult<Vec<u64>> {
    let mut m = modulus.to_vec();
    m = trim(m);
    if m.is_empty() {
        return Err(PyValueError::new_err("modulus polynomial is zero"));
    }
    let md = (m.len() as i64) - 1;
    let mlc = *m.last().unwrap();
    if mlc == 0 {
        return Err(PyValueError::new_err(
            "invalid modulus (leading coefficient 0)",
        ));
    }
    // We only support fields, so we require modulus to be monic for simplest reduction
    // TODO: generalization
    if mlc != 1 {
        return Err(PyValueError::new_err(
            "for now, modulus must be monic (leading coefficient 1)",
        ));
    }

    a = trim(a);
    while !a.is_empty() {
        let ad = (a.len() as i64) - 1;
        if ad < md {
            break;
        }
        let shift = (ad - md) as usize;
        let factor = *a.last().unwrap(); // since modulus is monic
        // a -= factor * x^shift * modulus
        for i in 0..m.len() {
            let idx = i + shift;
            let subtr = ((factor as u128 * m[i] as u128) % (p as u128)) as u64;
            let cur = a[idx];
            a[idx] = if cur >= subtr {
                cur - subtr
            } else {
                p - (subtr - cur)
            };
        }
        a = trim(a);
    }
    Ok(a)
}

#[pymethods]
impl Fq {
    /// Construct GF(p^k) as Fp[x]/(f). `modulus_coeffs` is low->high coeffs for f.
    /// For now, f must be monic and irreducible, which can be checked via PolyFp.is_irreducible().
    #[new]
    pub fn new(p: u64, modulus_coeffs: Vec<i128>) -> PyResult<Self> {
        let modulus = PolyFp::new(p, modulus_coeffs)?;
        if modulus.degree() < 1 {
            return Err(PyValueError::new_err("modulus degree must be >= 1"));
        }
        let mc = modulus.coeffs();
        if *mc.last().unwrap() != 1 {
            return Err(PyValueError::new_err(
                "for now, modulus must be monic (leading coefficient 1)",
            ));
        }
        Ok(Self {
            p,
            k: mc.len() - 1,
            modulus,
        })
    }

    /// Return the base characteristic.
    pub fn p(&self) -> u64 {
        self.p
    }
    /// Return the extension degree.
    pub fn degree(&self) -> usize {
        self.k
    }

    /// Return the field size `p^k`.
    pub fn size(&self) -> u64 {
        self.p.pow(self.k as u32)
    }

    /// Return the modulus polynomial coefficients.
    pub fn modulus_coeffs(&self) -> Vec<u64> {
        self.modulus.coeffs()
    }

    /// Create an element from polynomial coefficients and reduce modulo the field modulus.
    pub fn elem(&self, coeffs: Vec<i128>) -> PyResult<FqElem> {
        let p = self.p as i128;
        let mut c: Vec<u64> = coeffs
            .into_iter()
            .map(|x| {
                let mut r = x % p;
                if r < 0 {
                    r += p;
                }
                r as u64
            })
            .collect();
        c = trim(c);
        let reduced = poly_mod(self.p, c, &self.modulus.coeffs())?;
        Ok(FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: reduced,
        })
    }

    /// Return the additive identity.
    pub fn zero(&self) -> FqElem {
        FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: vec![],
        }
    }

    /// Return the multiplicative identity.
    pub fn one(&self) -> FqElem {
        FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: vec![1],
        }
    }

    /// Return `a + b`.
    pub fn add(&self, a: &FqElem, b: &FqElem) -> PyResult<FqElem> {
        self.check(a)?;
        self.check(b)?;
        Ok(FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: poly_add(self.p, &a.coeffs, &b.coeffs),
        })
    }

    /// Return `-a`.
    pub fn neg(&self, a: &FqElem) -> PyResult<FqElem> {
        self.check(a)?;
        Ok(FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: poly_neg(self.p, &a.coeffs),
        })
    }

    /// Return `a - b`.
    pub fn sub(&self, a: &FqElem, b: &FqElem) -> PyResult<FqElem> {
        self.add(a, &self.neg(b)?)
    }

    /// Return `a * b`, reduced modulo the field modulus.
    pub fn mul(&self, a: &FqElem, b: &FqElem) -> PyResult<FqElem> {
        self.check(a)?;
        self.check(b)?;
        let prod = poly_mul(self.p, &a.coeffs, &b.coeffs);
        let reduced = poly_mod(self.p, prod, &self.modulus.coeffs())?;
        Ok(FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: reduced,
        })
    }

    /// Return the permutation of field elements induced by `x -> x + b`.
    /// Uses the default enumeration bound `max_size = 4096`.
    pub fn add_perm(&self, b: &FqElem) -> PyResult<Perm> {
        self.add_perm_with_limit(b, 4096)
    }

    /// Return the permutation of field elements induced by `x -> x + b`.
    pub fn add_perm_with_limit(&self, b: &FqElem, max_size: u64) -> PyResult<Perm> {
        self.check(b)?;
        let elems = self.elements(Some(max_size))?;
        let mut images: Vec<usize> = Vec::with_capacity(elems.len());

        for x in &elems {
            let y = self.add(x, b)?;
            let idx = elems
                .iter()
                .position(|e| e.coeffs == y.coeffs)
                .ok_or_else(|| PyValueError::new_err("failed to index translated element"))?;
            images.push(idx);
        }

        Ok(Perm::new(elems.len(), images)?)
    }

    /// Return the permutation of field elements induced by `x -> a*x`.
    /// Uses the default enumeration bound `max_size = 4096`.
    pub fn mul_perm(&self, a: &FqElem) -> PyResult<Perm> {
        self.mul_perm_with_limit(a, 4096)
    }

    /// Return the permutation of field elements induced by `x -> a*x`.
    pub fn mul_perm_with_limit(&self, a: &FqElem, max_size: u64) -> PyResult<Perm> {
        self.check(a)?;
        if a.coeffs.is_empty() {
            return Err(PyValueError::new_err(
                "0 is not invertible, so x -> a*x is not a permutation",
            ));
        }
        let elems = self.elements(Some(max_size))?;
        let mut images: Vec<usize> = Vec::with_capacity(elems.len());

        for x in &elems {
            let y = self.mul(a, x)?;
            let idx = elems
                .iter()
                .position(|e| e.coeffs == y.coeffs)
                .ok_or_else(|| PyValueError::new_err("failed to index multiplied element"))?;
            images.push(idx);
        }

        Ok(Perm::new(elems.len(), images)?)
    }

    /// Return the affine permutation induced by `x -> a*x + b`.
    /// Uses the default enumeration bound `max_size = 4096`.
    pub fn affine_perm(&self, a: &FqElem, b: &FqElem) -> PyResult<Perm> {
        self.affine_perm_with_limit(a, b, 4096)
    }

    /// Return the affine permutation induced by `x -> a*x + b`.
    pub fn affine_perm_with_limit(&self, a: &FqElem, b: &FqElem, max_size: u64) -> PyResult<Perm> {
        self.check(a)?;
        self.check(b)?;
        if a.coeffs.is_empty() {
            return Err(PyValueError::new_err(
                "leading coefficient must be nonzero for x -> a*x + b to be a permutation",
            ));
        }
        let elems = self.elements(Some(max_size))?;
        let mut images: Vec<usize> = Vec::with_capacity(elems.len());

        for x in &elems {
            let ax = self.mul(a, x)?;
            let y = self.add(&ax, b)?;
            let idx = elems
                .iter()
                .position(|e| e.coeffs == y.coeffs)
                .ok_or_else(|| PyValueError::new_err("failed to index affine image"))?;
            images.push(idx);
        }

        Ok(Perm::new(elems.len(), images)?)
    }

    /// Return the multiplicative order of a nonzero element.
    pub fn mul_order(&self, a: &FqElem) -> PyResult<u64> {
        self.check(a)?;
        if a.coeffs.is_empty() {
            return Err(PyValueError::new_err(
                "0 is not in the multiplicative group",
            ));
        }

        let one = self.one();
        let mut acc = one.clone();
        let max = self.size() - 1; // q-1

        for k in 1..=max {
            acc = self.mul(&acc, a)?;
            if acc.coeffs == one.coeffs {
                return Ok(k);
            }
        }
        Err(PyValueError::new_err(
            "failed to find multiplicative order (unexpected)",
        ))
    }

    /// Return the multiplicative inverse of a nonzero element.
    pub fn inv(&self, a: &FqElem) -> PyResult<FqElem> {
        self.check(a)?;
        if a.coeffs.is_empty() {
            return Err(PyZeroDivisionError::new_err("0 has no inverse"));
        }
        let aa = PolyFp::new(self.p, a.coeffs.iter().map(|&x| x as i128).collect())?;
        let (g, s, _t) = aa.egcd(&self.modulus)?;
        if !g.is_one() {
            return Err(PyValueError::new_err(
                "element not invertible (unexpected if modulus irreducible)",
            ));
        }
        let reduced = poly_mod(self.p, s.coeffs(), &self.modulus.coeffs())?;
        Ok(FqElem {
            p: self.p,
            modulus_coeffs: self.modulus.coeffs(),
            coeffs: reduced,
        })
    }

    /// Return `a^exp`, allowing negative exponents for nonzero elements.
    pub fn pow(&self, a: &FqElem, exp: i128) -> PyResult<FqElem> {
        self.check(a)?;
        if exp == 0 {
            return Ok(self.one());
        }
        if exp < 0 {
            let inv = self.inv(a)?;
            return self.pow(&inv, -exp);
        }
        let mut e = exp as u64;
        let mut result = self.one();
        let mut base = a.clone();
        while e > 0 {
            if (e & 1) == 1 {
                result = self.mul(&result, &base)?;
            }
            e >>= 1;
            if e > 0 {
                base = self.mul(&base, &base)?;
            }
        }
        Ok(result)
    }

    /// Return up to `limit` primitive elements (generators of GF(q)^*).
    /// Uses defaults: max_size=4096, limit=64.
    pub fn primitive_elements(&self) -> PyResult<Vec<FqElem>> {
        self.primitive_elements_with_limit(4096, 64)
    }

    /// Return up to `limit` primitive elements, but refuse to enumerate if q > max_size.
    pub fn primitive_elements_with_limit(
        &self,
        max_size: u64,
        limit: usize,
    ) -> PyResult<Vec<FqElem>> {
        let q = self.size();
        if q > max_size {
            return Err(PyValueError::new_err(
                "field too large to enumerate (increase max_size)",
            ));
        }
        if q < 2 {
            return Ok(vec![]);
        }

        let target = q - 1;
        let elems = self.elements(Some(max_size))?; // includes 0
        let mut out: Vec<FqElem> = Vec::new();

        for a in elems.into_iter() {
            if a.coeffs.is_empty() {
                continue;
            } // skip 0
            let ord = self.mul_order(&a)?;
            if ord == target {
                out.push(a);
                if out.len() >= limit {
                    break;
                }
            }
        }
        Ok(out)
    }

    /// Enumerate all elements if size <= max_size
    pub fn elements(&self, max_size: Option<u64>) -> PyResult<Vec<FqElem>> {
        let max = max_size.unwrap_or(4096);
        let sz = self.size();
        if sz > max {
            return Err(PyValueError::new_err(
                "field too large to enumerate (increase max_size)",
            ));
        }
        let k = self.k;
        let p = self.p;

        let mut out = Vec::with_capacity(sz as usize);
        for mut t in 0..sz {
            let mut coeffs = vec![0u64; k];
            for i in 0..k {
                coeffs[i] = t % p;
                t /= p;
            }
            let coeffs = trim(coeffs);
            out.push(FqElem {
                p: self.p,
                modulus_coeffs: self.modulus.coeffs(),
                coeffs,
            });
        }
        Ok(out)
    }

    fn check(&self, a: &FqElem) -> PyResult<()> {
        if a.p != self.p || a.modulus_coeffs != self.modulus.coeffs() {
            Err(PyValueError::new_err("Mixed parents (different field)"))
        } else {
            Ok(())
        }
    }
}

#[pymethods]
impl FqElem {
    /// Return the reduced coefficient vector.
    pub fn coeffs(&self) -> Vec<u64> {
        self.coeffs.clone()
    }

    /// Return whether this element is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Return `False` only for zero.
    pub fn __bool__(&self) -> bool {
        !self.coeffs.is_empty()
    }

    /// Return a polynomial-style representation in `GF(p^k)`.
    pub fn __repr__(&self) -> String {
        let k = self.modulus_coeffs.len() - 1;
        if self.coeffs.is_empty() {
            return format!("0 in GF({}^{})", self.p, k);
        }

        // Build "1 + x^2 + x^5" style string
        let mut terms: Vec<String> = Vec::new();
        for (i, &c) in self.coeffs.iter().enumerate() {
            if c == 0 {
                continue;
            }
            let term = match (c, i) {
                (1, 0) => "1".to_string(),
                (_, 0) => format!("{}", c),
                (1, 1) => "x".to_string(),
                (_, 1) => format!("{}*x", c),
                (1, _) => format!("x^{}", i),
                (_, _) => format!("{}*x^{}", c, i),
            };
            terms.push(term);
        }
        if terms.is_empty() {
            return format!("0 in GF({}^{})", self.p, k);
        }
        format!("{} in GF({}^{})", terms.join(" + "), self.p, k)
    }

    /// Support `==` and `!=` for elements from the same field.
    pub fn __richcmp__(&self, other: &FqElem, op: pyo3::basic::CompareOp) -> PyResult<bool> {
        let eq = self.p == other.p
            && self.modulus_coeffs == other.modulus_coeffs
            && self.coeffs == other.coeffs;
        match op {
            pyo3::basic::CompareOp::Eq => Ok(eq),
            pyo3::basic::CompareOp::Ne => Ok(!eq),
            _ => Err(PyValueError::new_err("Only == and != supported")),
        }
    }

    // ---------- core helpers ----------
    fn check_same_parent(&self, other: &FqElem) -> PyResult<()> {
        if self.p != other.p || self.modulus_coeffs != other.modulus_coeffs {
            Err(PyValueError::new_err("Mixed parents (different field)"))
        } else {
            Ok(())
        }
    }

    fn canon_int(&self, k: i128) -> u64 {
        let m = self.p as i128;
        let mut r = k % m;
        if r < 0 {
            r += m;
        }
        r as u64
    }

    fn make_elem(&self, coeffs: Vec<u64>) -> FqElem {
        FqElem {
            p: self.p,
            modulus_coeffs: self.modulus_coeffs.clone(),
            coeffs: trim(coeffs),
        }
    }

    fn reduce(&self, coeffs: Vec<u64>) -> PyResult<FqElem> {
        let reduced = poly_mod(self.p, trim(coeffs), &self.modulus_coeffs)?;
        Ok(self.make_elem(reduced))
    }

    fn inv_internal(&self) -> PyResult<FqElem> {
        if self.coeffs.is_empty() {
            return Err(PyZeroDivisionError::new_err("0 has no inverse"));
        }
        let a = PolyFp::new(self.p, self.coeffs.iter().map(|&x| x as i128).collect())?;
        let f = PolyFp::new(
            self.p,
            self.modulus_coeffs.iter().map(|&x| x as i128).collect(),
        )?;
        let (g, s, _t) = a.egcd(&f)?;
        if !g.is_one() {
            return Err(PyValueError::new_err(
                "element not invertible (unexpected if modulus irreducible)",
            ));
        }
        self.reduce(s.coeffs())
    }

    fn pow_internal(&self, exp: i128) -> PyResult<FqElem> {
        if exp == 0 {
            return Ok(self.make_elem(vec![1]));
        }
        if exp < 0 {
            let inv = self.inv_internal()?;
            return inv.pow_internal(-exp);
        }
        let mut e = exp as u64;
        let mut result = self.make_elem(vec![1]);
        let mut base = self.clone();

        while e > 0 {
            if (e & 1) == 1 {
                result = result.__mul__(&base)?;
            }
            e >>= 1;
            if e > 0 {
                base = base.__mul__(&base)?;
            }
        }
        Ok(result)
    }

    /// Return the field trace down to the prime subfield.
    pub fn trace(&self) -> PyResult<FqElem> {
        // Tr(a) = sum_{i=0..k-1} a^(p^i)
        let k = self.modulus_coeffs.len() - 1;
        let mut acc = self.make_elem(vec![]); // 0
        for i in 0..k {
            let e = (self.p as i128).pow(i as u32); // p^i
            let term = self.pow_internal(e)?;
            acc = acc.__add__(&term)?;
        }
        Ok(acc)
    }

    /// Return the field norm down to the prime subfield.
    pub fn norm(&self) -> PyResult<FqElem> {
        // N(a) = product_{i=0..k-1} a^(p^i)
        let k = self.modulus_coeffs.len() - 1;
        let mut acc = self.make_elem(vec![1]); // 1
        for i in 0..k {
            let e = (self.p as i128).pow(i as u32); // p^i
            let term = self.pow_internal(e)?;
            acc = acc.__mul__(&term)?;
        }
        Ok(acc)
    }

    /// Return the multiplicative inverse of this nonzero element.
    pub fn inv(&self) -> PyResult<FqElem> {
        self.inv_internal()
    }

    /// Return the additive inverse.
    pub fn __neg__(&self) -> FqElem {
        let p = self.p;
        let coeffs = trim(
            self.coeffs
                .iter()
                .map(|&c| if c == 0 { 0 } else { p - c })
                .collect(),
        );
        self.make_elem(coeffs)
    }

    /// Add two extension-field elements.
    pub fn __add__(&self, other: &FqElem) -> PyResult<FqElem> {
        self.check_same_parent(other)?;
        Ok(self.make_elem(poly_add(self.p, &self.coeffs, &other.coeffs)))
    }

    /// Subtract two extension-field elements.
    pub fn __sub__(&self, other: &FqElem) -> PyResult<FqElem> {
        self.check_same_parent(other)?;
        Ok(self.make_elem(poly_add(
            self.p,
            &self.coeffs,
            &poly_neg(self.p, &other.coeffs),
        )))
    }

    /// Multiply two extension-field elements.
    pub fn __mul__(&self, other: &FqElem) -> PyResult<FqElem> {
        self.check_same_parent(other)?;
        let prod = poly_mul(self.p, &self.coeffs, &other.coeffs);
        self.reduce(prod)
    }

    /// Divide by another nonzero extension-field element.
    pub fn __truediv__(&self, other: &FqElem) -> PyResult<FqElem> {
        self.check_same_parent(other)?;
        let inv = other.inv_internal()?;
        self.__mul__(&inv)
    }

    /// Ternary pow protocol: pow(a, exp, mod). For `a ** exp`, mod is None.
    pub fn __pow__(&self, exp: i128, modulo: Option<Py<PyAny>>) -> PyResult<FqElem> {
        if modulo.is_some() {
            return Err(PyValueError::new_err(
                "3-argument pow(a, e, m) is not supported for FqElem",
            ));
        }
        self.pow_internal(exp)
    }

    /// Support `int + FqElem`.
    pub fn __radd__(&self, other: i128) -> PyResult<FqElem> {
        let c = self.canon_int(other);
        self.__add__(&self.make_elem(vec![c]))
    }

    /// Support `int - FqElem`.
    pub fn __rsub__(&self, other: i128) -> PyResult<FqElem> {
        let c = self.canon_int(other);
        self.make_elem(vec![c]).__sub__(self)
    }

    /// Support `int * FqElem`.
    pub fn __rmul__(&self, other: i128) -> PyResult<FqElem> {
        let c = self.canon_int(other);
        self.__mul__(&self.make_elem(vec![c]))
    }

    /// Support `int / FqElem`.
    pub fn __rtruediv__(&self, other: i128) -> PyResult<FqElem> {
        let c = self.canon_int(other);
        self.make_elem(vec![c]).__truediv__(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fq_inverse_via_egcd_multiplies_to_one() {
        let k = Fq::new(2, vec![1, 1, 0, 1]).unwrap();
        let a = k.elem(vec![1, 0, 1]).unwrap();
        let inv = k.inv(&a).unwrap();
        let prod = k.mul(&a, &inv).unwrap();
        assert_eq!(prod.coeffs(), vec![1]);
    }

    #[test]
    fn fq_negative_power_uses_inverse_path() {
        let k = Fq::new(2, vec![1, 1, 0, 1]).unwrap();
        let a = k.elem(vec![0, 1]).unwrap();
        let inv = k.inv(&a).unwrap();
        let pow = k.pow(&a, -1).unwrap();
        assert_eq!(pow.coeffs(), inv.coeffs());
    }

    #[test]
    fn fq_add_mul_and_affine_perms_exist() {
        let k = Fq::new(2, vec![1, 1, 0, 1]).unwrap();
        let a = k.elem(vec![0, 1]).unwrap();
        let b = k.elem(vec![1]).unwrap();

        let add = k.add_perm(&b).unwrap();
        let mul = k.mul_perm(&a).unwrap();
        let affine = k.affine_perm(&a, &b).unwrap();

        assert_eq!(add.n(), 8);
        assert_eq!(mul.n(), 8);
        assert_eq!(add.compose(&mul).unwrap(), affine);
    }

    #[test]
    fn fq_mul_perm_rejects_zero() {
        let k = Fq::new(2, vec![1, 1, 0, 1]).unwrap();
        assert!(k.mul_perm(&k.zero()).is_err());
        assert!(k.affine_perm(&k.zero(), &k.one()).is_err());
    }
}

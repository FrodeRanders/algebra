use pyo3::prelude::*;
use pyo3::exceptions::{PyValueError, PyZeroDivisionError};

use crate::arith::egcd::inv_mod_i128;
use crate::arith::prime::is_prime_u64;

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct Fp {
    p: u64,
}

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct FpElem {
    p: u64,
    v: u64,
}

fn canon_mod_i128(k: i128, p: u64) -> u64 {
    let m = p as i128;
    let mut r = k % m;
    if r < 0 { r += m; }
    r as u64
}

#[pymethods]
impl Fp {
    /// Create a prime field GF(p). Requires p to be prime.
    #[new]
    pub fn new(p: u64) -> PyResult<Self> {
        if !is_prime_u64(p) {
            return Err(PyValueError::new_err("p must be prime for Fp(p)"));
        }
        Ok(Self { p })
    }

    pub fn modulus(&self) -> u64 { self.p }

    pub fn zero(&self) -> FpElem { FpElem { p: self.p, v: 0 } }
    pub fn one(&self) -> FpElem { FpElem { p: self.p, v: 1 % self.p } }

    /// Create an element from (possibly negative) integer.
    pub fn elem(&self, k: i128) -> FpElem {
        FpElem { p: self.p, v: canon_mod_i128(k, self.p) }
    }

    /// Enumerate all elements (0..p-1).
    pub fn elements(&self) -> Vec<FpElem> {
        (0..self.p).map(|v| FpElem { p: self.p, v }).collect()
    }

    /// Enumerate all nonzero elements (1..p-1).
    pub fn nonzero_elements(&self) -> Vec<FpElem> {
        (1..self.p).map(|v| FpElem { p: self.p, v }).collect()
    }

    pub fn add(&self, a: &FpElem, b: &FpElem) -> PyResult<FpElem> {
        self.check(a)?;
        self.check(b)?;
        Ok(FpElem { p: self.p, v: (a.v + b.v) % self.p })
    }

    pub fn neg(&self, a: &FpElem) -> PyResult<FpElem> {
        self.check(a)?;
        Ok(if a.v == 0 {
            a.clone()
        } else {
            FpElem { p: self.p, v: self.p - a.v }
        })
    }

    pub fn sub(&self, a: &FpElem, b: &FpElem) -> PyResult<FpElem> {
        self.add(a, &self.neg(b)?)
    }

    pub fn mul(&self, a: &FpElem, b: &FpElem) -> PyResult<FpElem> {
        self.check(a)?;
        self.check(b)?;
        let prod = (a.v as u128) * (b.v as u128);
        Ok(FpElem { p: self.p, v: (prod % (self.p as u128)) as u64 })
    }

    /// Multiplicative inverse in GF(p). Raises ZeroDivisionError for 0.
    pub fn inv(&self, a: &FpElem) -> PyResult<FpElem> {
        self.check(a)?;
        if a.v == 0 {
            return Err(PyZeroDivisionError::new_err("0 has no multiplicative inverse in a field"));
        }
        let inv = inv_mod_i128(a.v as i128, self.p as i128)
            .expect("inv_mod should exist in prime field");
        Ok(FpElem { p: self.p, v: inv as u64 })
    }

    /// Exponentiation a^e (supports negative exponents for nonzero a).
    pub fn pow(&self, a: &FpElem, e: i128) -> PyResult<FpElem> {
        self.check(a)?;
        if e == 0 {
            return Ok(self.one());
        }
        if e < 0 {
            if a.v == 0 {
                return Err(PyZeroDivisionError::new_err("0 cannot be raised to a negative power"));
            }
            let inv_a = self.inv(a)?;
            return self.pow(&inv_a, -e);
        }
        Ok(pow_u64_mod(a.v, e as u64, self.p))
    }

    /// Multiplicative order of a in GF(p)^* (brute force; p should be small-ish).
    /// Raises ValueError if a = 0.
    pub fn mul_order(&self, a: &FpElem) -> PyResult<u64> {
        self.check(a)?;
        if a.v == 0 {
            return Err(PyValueError::new_err("0 is not in the multiplicative group"));
        }
        // group size is p-1
        let one = self.one();
        let mut acc = one.clone();
        for k in 1..=self.p - 1 {
            acc = self.mul(&acc, a)?;
            if acc.v == one.v {
                return Ok(k);
            }
        }
        // Should never happen in a field for nonzero element.
        Err(PyValueError::new_err("failed to find multiplicative order (unexpected)"))
    }

    /// Sanity utility: multiplicative structure is abelian (always true for fields).
    pub fn is_abelian_mul(&self) -> bool { true }

    fn check(&self, a: &FpElem) -> PyResult<()> {
        if a.p != self.p {
            Err(PyValueError::new_err("Mixed moduli / parents"))
        } else {
            Ok(())
        }
    }

    /// Euler criterion: a is a quadratic residue iff a^((p-1)/2) == 1 (and a != 0).
    /// Returns True for a=0 by convention.
    pub fn is_quadratic_residue(&self, a: &FpElem) -> PyResult<bool> {
        self.check(a)?;
        if a.v == 0 {
            return Ok(true);
        }
        let e = ((self.p - 1) / 2) as i128;
        let t = self.pow(a, e)?;
        Ok(t.v == 1)
    }

    /// Brute-force square root in Fp (good for small primes).
    /// Returns one square root if it exists, else None.
    pub fn sqrt(&self, a: &FpElem) -> PyResult<Option<FpElem>> {
        self.check(a)?;
        // For course primes, brute force is fine.
        for x in 0..self.p {
            let xe = FpElem { p: self.p, v: x };
            let sq = self.mul(&xe, &xe)?;
            if sq.v == a.v {
                return Ok(Some(xe));
            }
        }
        Ok(None)
    }

    /// Discrete log in Fp* (brute force): find x s.t. g^x = h.
    /// Only for small p. Raises ValueError if g=0 or h=0 or no solution.
    pub fn dlog(&self, g: &FpElem, h: &FpElem) -> PyResult<u64> {
        self.check(g)?;
        self.check(h)?;
        if g.v == 0 || h.v == 0 {
            return Err(PyValueError::new_err("dlog is defined in the multiplicative group (nonzero elements)"));
        }
        let mut acc = self.one();
        for x in 0..=(self.p - 2) {
            if acc.v == h.v {
                return Ok(x);
            }
            acc = self.mul(&acc, g)?;
        }
        Err(PyValueError::new_err("No discrete log found (g may not generate the subgroup containing h)"))
    }
}

#[pymethods]
impl FpElem {
    pub fn value(&self) -> u64 { self.v }

    pub fn __repr__(&self) -> String {
        format!("{} (mod {})", self.v, self.p)
    }

    pub fn __int__(&self) -> u64 { self.v }

    pub fn __richcmp__(&self, other: &FpElem, op: pyo3::basic::CompareOp) -> PyResult<bool> {
        if self.p != other.p {
            return Ok(false);
        }
        match op {
            pyo3::basic::CompareOp::Eq => Ok(self.v == other.v),
            pyo3::basic::CompareOp::Ne => Ok(self.v != other.v),
            _ => Err(PyValueError::new_err("Only == and != supported")),
        }
    }

    // ---------- element-level helpers ----------
    fn check_same_p(&self, other: &FpElem) -> PyResult<()> {
        if self.p != other.p {
            Err(PyValueError::new_err("Mixed moduli / parents"))
        } else {
            Ok(())
        }
    }

    fn canon_int(&self, k: i128) -> u64 {
        let m = self.p as i128;
        let mut r = k % m;
        if r < 0 { r += m; }
        r as u64
    }

    fn add_v(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }

    fn mul_v(&self, a: u64, b: u64) -> u64 {
        let prod = (a as u128) * (b as u128);
        (prod % (self.p as u128)) as u64
    }

    fn inv_v(&self, a: u64) -> PyResult<u64> {
        if a == 0 {
            return Err(PyZeroDivisionError::new_err("0 has no multiplicative inverse in a field"));
        }
        let inv = inv_mod_i128(a as i128, self.p as i128)
            .expect("inverse should exist for nonzero element in prime field");
        Ok(inv as u64)
    }

    fn pow_v(&self, base: u64, exp: i128) -> PyResult<u64> {
        if exp == 0 {
            return Ok(1 % self.p);
        }
        if exp < 0 {
            if base == 0 {
                return Err(PyZeroDivisionError::new_err("0 cannot be raised to a negative power"));
            }
            let inv = self.inv_v(base)?;
            return self.pow_v(inv, -exp);
        }
        let mut e = exp as u64;
        let mut result: u64 = 1 % self.p;
        let mut b: u64 = base % self.p;

        while e > 0 {
            if (e & 1) == 1 {
                result = self.mul_v(result, b);
            }
            e >>= 1;
            if e > 0 {
                b = self.mul_v(b, b);
            }
        }
        Ok(result)
    }

    pub fn mul_order(&self) -> PyResult<u64> {
        if self.v == 0 {
            return Err(PyValueError::new_err("0 is not in the multiplicative group"));
        }
        let one = 1 % self.p;
        let mut acc: u64 = one;

        for k in 1..=self.p - 1 {
            let prod = (acc as u128) * (self.v as u128);
            acc = (prod % (self.p as u128)) as u64;
            if acc == one {
                return Ok(k);
            }
        }
        Err(PyValueError::new_err("failed to find multiplicative order (unexpected)"))
    }

    // ---------- nice API ----------
    pub fn inv(&self) -> PyResult<FpElem> {
        Ok(FpElem { p: self.p, v: self.inv_v(self.v)? })
    }

    // ---------- operator overloading ----------
    pub fn __neg__(&self) -> FpElem {
        if self.v == 0 { self.clone() }
        else { FpElem { p: self.p, v: self.p - self.v } }
    }

    pub fn __add__(&self, other: &FpElem) -> PyResult<FpElem> {
        self.check_same_p(other)?;
        Ok(FpElem { p: self.p, v: self.add_v(self.v, other.v) })
    }

    pub fn __sub__(&self, other: &FpElem) -> PyResult<FpElem> {
        self.check_same_p(other)?;
        let neg = if other.v == 0 { 0 } else { self.p - other.v };
        Ok(FpElem { p: self.p, v: self.add_v(self.v, neg) })
    }

    pub fn __mul__(&self, other: &FpElem) -> PyResult<FpElem> {
        self.check_same_p(other)?;
        Ok(FpElem { p: self.p, v: self.mul_v(self.v, other.v) })
    }

    pub fn __truediv__(&self, other: &FpElem) -> PyResult<FpElem> {
        self.check_same_p(other)?;
        let inv = self.inv_v(other.v)?;
        Ok(FpElem { p: self.p, v: self.mul_v(self.v, inv) })
    }

    /// Supports `a ** exp` and `pow(a, exp)`. Rejects `pow(a, exp, mod)` if mod is not None.
    pub fn __pow__(&self, exp: i128, modulo: Option<Py<PyAny>>) -> PyResult<FpElem> {        if modulo.is_some() {
            return Err(PyValueError::new_err(
                "3-argument pow(a, e, m) is not supported for FpElem",
            ));
        }
        Ok(FpElem { p: self.p, v: self.pow_v(self.v, exp)? })
    }

    // ---------- mixed int ops (int OP elem) ----------
    pub fn __radd__(&self, other: i128) -> PyResult<FpElem> {
        let ov = self.canon_int(other);
        Ok(FpElem { p: self.p, v: self.add_v(ov, self.v) })
    }

    pub fn __rsub__(&self, other: i128) -> PyResult<FpElem> {
        let ov = self.canon_int(other);
        let neg_self = if self.v == 0 { 0 } else { self.p - self.v };
        Ok(FpElem { p: self.p, v: self.add_v(ov, neg_self) })
    }

    pub fn __rmul__(&self, other: i128) -> PyResult<FpElem> {
        let ov = self.canon_int(other);
        Ok(FpElem { p: self.p, v: self.mul_v(ov, self.v) })
    }

    pub fn __rtruediv__(&self, other: i128) -> PyResult<FpElem> {
        let ov = self.canon_int(other);
        let inv = self.inv_v(self.v)?;
        Ok(FpElem { p: self.p, v: self.mul_v(ov, inv) })
    }
}

fn pow_u64_mod(base: u64, mut exp: u64, p: u64) -> FpElem {
    let mut result: u64 = 1 % p;
    let mut b: u64 = base % p;

    while exp > 0 {
        if (exp & 1) == 1 {
            let prod = (result as u128) * (b as u128);
            result = (prod % (p as u128)) as u64;
        }
        exp >>= 1;
        if exp > 0 {
            let sq = (b as u128) * (b as u128);
            b = (sq % (p as u128)) as u64;
        }
    }
    FpElem { p, v: result }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fp_requires_prime() {
        assert!(Fp::new(7).is_ok());
        assert!(Fp::new(12).is_err());
        assert!(Fp::new(1).is_err());
    }

    #[test]
    fn fp_canonicalizes_negatives() {
        let f = Fp::new(7).unwrap();
        assert_eq!(f.elem(-1).v, 6);
        assert_eq!(f.elem(-8).v, 6);
        assert_eq!(f.elem(8).v, 1);
    }

    #[test]
    fn fp_inv_works() {
        let f = Fp::new(7).unwrap();
        let a = f.elem(3);
        let inv = f.inv(&a).unwrap();
        let prod = f.mul(&a, &inv).unwrap();
        assert_eq!(prod.v, 1);
    }

    #[test]
    fn fp_pow_negative() {
        let f = Fp::new(11).unwrap();
        let a = f.elem(2);
        // 2^-1 mod 11 = 6, so 2^-2 = 6^2 = 36 mod 11 = 3
        let x = f.pow(&a, -2).unwrap();
        assert_eq!(x.v, 3);
    }

    #[test]
    fn fp_mul_order() {
        let f = Fp::new(7).unwrap();
        let a = f.elem(3);
        // In F7*, 3 has order 6 (primitive element)
        assert_eq!(f.mul_order(&a).unwrap(), 6);
    }
}

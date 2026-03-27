use pyo3::exceptions::{PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::{Py, PyAny};

use crate::arith::egcd::inv_mod_i128;

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct Zn {
    n: u64,
}

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ZnElem {
    n: u64,
    v: u64,
}

fn canon_mod_i128(k: i128, n: u64) -> u64 {
    let m = n as i128;
    let mut r = k % m;
    if r < 0 {
        r += m;
    }
    r as u64
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

#[pymethods]
impl Zn {
    #[new]
    pub fn new(n: u64) -> PyResult<Self> {
        if n < 2 {
            return Err(PyValueError::new_err("n must be >= 2 for Zn(n)"));
        }
        Ok(Self { n })
    }

    pub fn modulus(&self) -> u64 {
        self.n
    }

    pub fn zero(&self) -> ZnElem {
        ZnElem { n: self.n, v: 0 }
    }

    pub fn one(&self) -> ZnElem {
        ZnElem {
            n: self.n,
            v: 1 % self.n,
        }
    }

    pub fn elem(&self, k: i128) -> ZnElem {
        ZnElem {
            n: self.n,
            v: canon_mod_i128(k, self.n),
        }
    }

    pub fn elements(&self) -> Vec<ZnElem> {
        (0..self.n).map(|v| ZnElem { n: self.n, v }).collect()
    }

    pub fn units(&self) -> Vec<ZnElem> {
        (0..self.n)
            .filter(|&v| gcd_u64(v, self.n) == 1)
            .map(|v| ZnElem { n: self.n, v })
            .collect()
    }

    pub fn is_unit(&self, a: &ZnElem) -> PyResult<bool> {
        self.check(a)?;
        Ok(gcd_u64(a.v, self.n) == 1)
    }

    pub fn add(&self, a: &ZnElem, b: &ZnElem) -> PyResult<ZnElem> {
        self.check(a)?;
        self.check(b)?;
        Ok(ZnElem {
            n: self.n,
            v: (a.v + b.v) % self.n,
        })
    }

    pub fn neg(&self, a: &ZnElem) -> PyResult<ZnElem> {
        self.check(a)?;
        Ok(if a.v == 0 {
            a.clone()
        } else {
            ZnElem {
                n: self.n,
                v: self.n - a.v,
            }
        })
    }

    pub fn sub(&self, a: &ZnElem, b: &ZnElem) -> PyResult<ZnElem> {
        self.add(a, &self.neg(b)?)
    }

    pub fn mul(&self, a: &ZnElem, b: &ZnElem) -> PyResult<ZnElem> {
        self.check(a)?;
        self.check(b)?;
        let prod = (a.v as u128) * (b.v as u128);
        Ok(ZnElem {
            n: self.n,
            v: (prod % (self.n as u128)) as u64,
        })
    }

    pub fn inv(&self, a: &ZnElem) -> PyResult<ZnElem> {
        self.check(a)?;
        if a.v == 0 {
            return Err(PyZeroDivisionError::new_err(
                "0 has no multiplicative inverse in Zn",
            ));
        }
        let inv = inv_mod_i128(a.v as i128, self.n as i128)
            .ok_or_else(|| PyValueError::new_err("element is not a unit in Zn"))?;
        Ok(ZnElem {
            n: self.n,
            v: inv as u64,
        })
    }

    pub fn pow(&self, a: &ZnElem, e: i128) -> PyResult<ZnElem> {
        self.check(a)?;
        if e == 0 {
            return Ok(self.one());
        }
        if e < 0 {
            let inv = self.inv(a)?;
            return self.pow(&inv, -e);
        }

        let mut exp = e as u64;
        let mut result: u64 = 1 % self.n;
        let mut base: u64 = a.v;

        while exp > 0 {
            if (exp & 1) == 1 {
                result = ((result as u128 * base as u128) % (self.n as u128)) as u64;
            }
            exp >>= 1;
            if exp > 0 {
                base = ((base as u128 * base as u128) % (self.n as u128)) as u64;
            }
        }

        Ok(ZnElem {
            n: self.n,
            v: result,
        })
    }

    fn check(&self, a: &ZnElem) -> PyResult<()> {
        if a.n != self.n {
            Err(PyValueError::new_err("Mixed moduli / parents"))
        } else {
            Ok(())
        }
    }
}

#[pymethods]
impl ZnElem {
    pub fn value(&self) -> u64 {
        self.v
    }

    pub fn is_zero(&self) -> bool {
        self.v == 0
    }

    pub fn is_unit(&self) -> bool {
        gcd_u64(self.v, self.n) == 1
    }

    pub fn __bool__(&self) -> bool {
        self.v != 0
    }

    pub fn __repr__(&self) -> String {
        format!("{} (mod {})", self.v, self.n)
    }

    pub fn __int__(&self) -> u64 {
        self.v
    }

    pub fn __richcmp__(&self, other: &ZnElem, op: pyo3::basic::CompareOp) -> PyResult<bool> {
        if self.n != other.n {
            return Ok(false);
        }
        match op {
            pyo3::basic::CompareOp::Eq => Ok(self.v == other.v),
            pyo3::basic::CompareOp::Ne => Ok(self.v != other.v),
            _ => Err(PyValueError::new_err("Only == and != supported")),
        }
    }

    fn check_same_n(&self, other: &ZnElem) -> PyResult<()> {
        if self.n != other.n {
            Err(PyValueError::new_err("Mixed moduli / parents"))
        } else {
            Ok(())
        }
    }

    fn canon_int(&self, k: i128) -> u64 {
        canon_mod_i128(k, self.n)
    }

    fn add_v(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.n
    }

    fn mul_v(&self, a: u64, b: u64) -> u64 {
        let prod = (a as u128) * (b as u128);
        (prod % (self.n as u128)) as u64
    }

    fn inv_v(&self, a: u64) -> PyResult<u64> {
        if a == 0 {
            return Err(PyZeroDivisionError::new_err(
                "0 has no multiplicative inverse in Zn",
            ));
        }
        let inv = inv_mod_i128(a as i128, self.n as i128)
            .ok_or_else(|| PyValueError::new_err("element is not a unit in Zn"))?;
        Ok(inv as u64)
    }

    fn pow_v(&self, base: u64, exp: i128) -> PyResult<u64> {
        if exp == 0 {
            return Ok(1 % self.n);
        }
        if exp < 0 {
            let inv = self.inv_v(base)?;
            return self.pow_v(inv, -exp);
        }

        let mut e = exp as u64;
        let mut result: u64 = 1 % self.n;
        let mut b: u64 = base % self.n;

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

    pub fn inv(&self) -> PyResult<ZnElem> {
        Ok(ZnElem {
            n: self.n,
            v: self.inv_v(self.v)?,
        })
    }

    pub fn __neg__(&self) -> ZnElem {
        if self.v == 0 {
            self.clone()
        } else {
            ZnElem {
                n: self.n,
                v: self.n - self.v,
            }
        }
    }

    pub fn __add__(&self, other: &ZnElem) -> PyResult<ZnElem> {
        self.check_same_n(other)?;
        Ok(ZnElem {
            n: self.n,
            v: self.add_v(self.v, other.v),
        })
    }

    pub fn __sub__(&self, other: &ZnElem) -> PyResult<ZnElem> {
        self.check_same_n(other)?;
        let neg = if other.v == 0 { 0 } else { self.n - other.v };
        Ok(ZnElem {
            n: self.n,
            v: self.add_v(self.v, neg),
        })
    }

    pub fn __mul__(&self, other: &ZnElem) -> PyResult<ZnElem> {
        self.check_same_n(other)?;
        Ok(ZnElem {
            n: self.n,
            v: self.mul_v(self.v, other.v),
        })
    }

    pub fn __truediv__(&self, other: &ZnElem) -> PyResult<ZnElem> {
        self.check_same_n(other)?;
        let inv = self.inv_v(other.v)?;
        Ok(ZnElem {
            n: self.n,
            v: self.mul_v(self.v, inv),
        })
    }

    pub fn __pow__(&self, exp: i128, modulo: Option<Py<PyAny>>) -> PyResult<ZnElem> {
        if modulo.is_some() {
            return Err(PyValueError::new_err(
                "3-argument pow(a, e, m) is not supported for ZnElem",
            ));
        }
        Ok(ZnElem {
            n: self.n,
            v: self.pow_v(self.v, exp)?,
        })
    }

    pub fn __radd__(&self, other: i128) -> PyResult<ZnElem> {
        let ov = self.canon_int(other);
        Ok(ZnElem {
            n: self.n,
            v: self.add_v(ov, self.v),
        })
    }

    pub fn __rsub__(&self, other: i128) -> PyResult<ZnElem> {
        let ov = self.canon_int(other);
        let neg_self = if self.v == 0 { 0 } else { self.n - self.v };
        Ok(ZnElem {
            n: self.n,
            v: self.add_v(ov, neg_self),
        })
    }

    pub fn __rmul__(&self, other: i128) -> PyResult<ZnElem> {
        let ov = self.canon_int(other);
        Ok(ZnElem {
            n: self.n,
            v: self.mul_v(ov, self.v),
        })
    }

    pub fn __rtruediv__(&self, other: i128) -> PyResult<ZnElem> {
        let ov = self.canon_int(other);
        let inv = self.inv_v(self.v)?;
        Ok(ZnElem {
            n: self.n,
            v: self.mul_v(ov, inv),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zn_requires_modulus_at_least_two() {
        assert!(Zn::new(2).is_ok());
        assert!(Zn::new(1).is_err());
    }

    #[test]
    fn zn_canonicalizes_negatives() {
        let z = Zn::new(8).unwrap();
        assert_eq!(z.elem(-1).v, 7);
        assert_eq!(z.elem(-9).v, 7);
        assert_eq!(z.elem(9).v, 1);
    }

    #[test]
    fn zn_units_are_detected() {
        let z = Zn::new(8).unwrap();
        let units: Vec<u64> = z.units().into_iter().map(|a| a.v).collect();
        assert_eq!(units, vec![1, 3, 5, 7]);
        assert!(z.is_unit(&z.elem(3)).unwrap());
        assert!(!z.is_unit(&z.elem(6)).unwrap());
    }

    #[test]
    fn zn_inverse_exists_only_for_units() {
        let z = Zn::new(10).unwrap();
        let a = z.elem(3);
        let inv = z.inv(&a).unwrap();
        assert_eq!(z.mul(&a, &inv).unwrap().v, 1);
        assert!(z.inv(&z.elem(2)).is_err());
    }

    #[test]
    fn zn_negative_power_requires_unit() {
        let z = Zn::new(9).unwrap();
        assert_eq!(z.pow(&z.elem(2), -1).unwrap().v, 5);
        assert!(z.pow(&z.elem(3), -1).is_err());
    }

    #[test]
    fn zn_elem_operator_mixing_works() {
        let z = Zn::new(12).unwrap();
        let a = z.elem(7);
        let b = z.elem(11);
        assert_eq!((a.__add__(&b).unwrap()).v, 6);
        assert_eq!((a.__mul__(&b).unwrap()).v, 5);
        assert_eq!(a.__radd__(5).unwrap().v, 0);
        assert_eq!(a.__rmul__(5).unwrap().v, 11);
    }
}

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::arith::egcd::inv_mod_i128;

/// A polynomial over `F_p` stored from low degree to high degree.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PolyFp {
    p: u64,
    coeffs: Vec<u64>, // low -> high
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

fn canon_coeffs(p: u64, coeffs: &[i128]) -> Vec<u64> {
    let m = p as i128;
    trim(
        coeffs
            .iter()
            .map(|&c| {
                let mut r = c % m;
                if r < 0 {
                    r += m;
                }
                r as u64
            })
            .collect(),
    )
}

fn inv_mod_p(a: u64, p: u64) -> PyResult<u64> {
    if a == 0 {
        return Err(PyValueError::new_err(
            "leading coefficient not invertible (0)",
        ));
    }
    let inv = inv_mod_i128(a as i128, p as i128)
        .ok_or_else(|| PyValueError::new_err("coefficient not invertible mod p"))?;
    Ok(inv as u64)
}

#[pymethods]
impl PolyFp {
    #[new]
    /// Construct a polynomial over `F_p`.
    pub fn new(p: u64, coeffs: Vec<i128>) -> PyResult<Self> {
        if p < 2 {
            return Err(PyValueError::new_err("p must be >= 2"));
        }
        Ok(Self {
            p,
            coeffs: canon_coeffs(p, &coeffs),
        })
    }

    /// Return the base characteristic.
    pub fn p(&self) -> u64 {
        self.p
    }

    /// Return the coefficient list.
    pub fn coeffs(&self) -> Vec<u64> {
        self.coeffs.clone()
    }

    /// Return the degree, or `-1` for the zero polynomial.
    pub fn degree(&self) -> i64 {
        if self.coeffs.is_empty() {
            -1
        } else {
            (self.coeffs.len() as i64) - 1
        }
    }

    /// Return whether the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Return whether the polynomial is the constant polynomial `1`.
    pub fn is_one(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0] == 1
    }

    /// Return the monic normalization of this polynomial.
    pub fn monic(&self) -> PyResult<PolyFp> {
        if self.is_zero() {
            return Err(PyValueError::new_err(
                "zero polynomial cannot be made monic",
            ));
        }
        let lc = *self.coeffs.last().unwrap();
        let inv = inv_mod_p(lc, self.p)?;
        Ok(self.scale(inv))
    }

    /// Return the sum of two polynomials.
    pub fn add(&self, other: &PolyFp) -> PyResult<PolyFp> {
        self.check(other)?;
        Ok(PolyFp {
            p: self.p,
            coeffs: trim(add_coeffs(self.p, &self.coeffs, &other.coeffs)),
        })
    }

    /// Return the additive inverse.
    pub fn neg(&self) -> PolyFp {
        let p = self.p;
        let coeffs = self
            .coeffs
            .iter()
            .map(|&c| if c == 0 { 0 } else { p - c })
            .collect();
        PolyFp {
            p,
            coeffs: trim(coeffs),
        }
    }

    /// Return the difference of two polynomials.
    pub fn sub(&self, other: &PolyFp) -> PyResult<PolyFp> {
        self.add(&other.neg())
    }

    /// Return the product of two polynomials.
    pub fn mul(&self, other: &PolyFp) -> PyResult<PolyFp> {
        self.check(other)?;
        Ok(PolyFp {
            p: self.p,
            coeffs: trim(mul_coeffs(self.p, &self.coeffs, &other.coeffs)),
        })
    }

    /// Multiply by a scalar modulo `p`.
    pub fn scale(&self, k: u64) -> PolyFp {
        let p = self.p;
        if k % p == 0 || self.is_zero() {
            return PolyFp { p, coeffs: vec![] };
        }
        let coeffs = self
            .coeffs
            .iter()
            .map(|&c| ((c as u128 * k as u128) % (p as u128)) as u64)
            .collect();
        PolyFp {
            p,
            coeffs: trim(coeffs),
        }
    }

    /// Long division: self = q*modulus + r
    pub fn div_rem(&self, modulus: &PolyFp) -> PyResult<(PolyFp, PolyFp)> {
        self.check(modulus)?;
        if modulus.is_zero() {
            return Err(PyValueError::new_err("division by zero polynomial"));
        }
        let p = self.p;
        let mut r = self.coeffs.clone();
        let mut q: Vec<u64> = vec![];

        let md = modulus.degree();
        let mut rdeg = if r.is_empty() {
            -1
        } else {
            (r.len() as i64) - 1
        };
        let mlc = *modulus.coeffs.last().unwrap();
        let mlc_inv = inv_mod_p(mlc, p)?;

        if rdeg < md {
            return Ok((PolyFp { p, coeffs: vec![] }, PolyFp { p, coeffs: trim(r) }));
        }

        q.resize((rdeg - md + 1) as usize, 0);

        while rdeg >= md && !r.is_empty() {
            let shift = (rdeg - md) as usize;
            let rlc = *r.last().unwrap();
            let factor = ((rlc as u128 * mlc_inv as u128) % (p as u128)) as u64;
            q[shift] = factor;

            // r -= factor * x^shift * modulus
            for i in 0..modulus.coeffs.len() {
                let idx = i + shift;
                let subtr = ((factor as u128 * modulus.coeffs[i] as u128) % (p as u128)) as u64;
                let cur = r[idx];
                r[idx] = if cur >= subtr {
                    cur - subtr
                } else {
                    p - (subtr - cur)
                };
            }
            r = trim(r);
            rdeg = if r.is_empty() {
                -1
            } else {
                (r.len() as i64) - 1
            };
        }

        Ok((PolyFp { p, coeffs: trim(q) }, PolyFp { p, coeffs: trim(r) }))
    }

    /// Return the remainder after division by `modulus`.
    pub fn modulo(&self, modulus: &PolyFp) -> PyResult<PolyFp> {
        let (_q, r) = self.div_rem(modulus)?;
        Ok(r)
    }

    /// Return the monic greatest common divisor.
    pub fn gcd(&self, other: &PolyFp) -> PyResult<PolyFp> {
        self.check(other)?;
        let mut a = self.clone();
        let mut b = other.clone();
        while !b.is_zero() {
            let r = a.modulo(&b)?;
            a = b;
            b = r;
        }
        if a.is_zero() { Ok(a) } else { a.monic() }
    }

    /// Return `(g, s, t)` such that `g = s*self + t*other`.
    pub fn egcd(&self, other: &PolyFp) -> PyResult<(PolyFp, PolyFp, PolyFp)> {
        self.check(other)?;
        let mut r0 = self.clone();
        let mut r1 = other.clone();
        let mut s0 = PolyFp::new(self.p, vec![1])?;
        let mut s1 = PolyFp::new(self.p, vec![0])?;
        let mut t0 = PolyFp::new(self.p, vec![0])?;
        let mut t1 = PolyFp::new(self.p, vec![1])?;

        while !r1.is_zero() {
            let (q, r) = r0.div_rem(&r1)?;
            let s = s0.sub(&q.mul(&s1)?)?;
            let t = t0.sub(&q.mul(&t1)?)?;
            r0 = r1;
            r1 = r;
            s0 = s1;
            s1 = s;
            t0 = t1;
            t1 = t;
        }
        Ok((r0, s0, t0))
    }

    /// Brute-force irreducibility check for small degrees.
    /// Tries all monic polynomials g with 1 <= deg(g) <= floor(deg(f)/2) and tests divisibility.
    pub fn is_irreducible(&self) -> PyResult<bool> {
        if self.is_zero() {
            return Ok(false);
        }
        let d = self.degree();
        if d <= 0 {
            return Ok(true); // constant or linear nonzero
        }
        // Ensure monic for a standard notion
        let f = self.monic()?;
        let max_g = d / 2;

        for deg_g in 1..=max_g {
            // enumerate all monic polynomials of degree deg_g over Fp
            let count = self.p.pow(deg_g as u32); // number of choices for lower coeffs
            for t in 0..count {
                let mut coeffs = vec![0u64; (deg_g as usize) + 1];
                coeffs[deg_g as usize] = 1;
                let mut x = t;
                for i in 0..(deg_g as usize) {
                    coeffs[i] = x % self.p;
                    x /= self.p;
                }
                let g = PolyFp {
                    p: self.p,
                    coeffs: trim(coeffs),
                };
                let (_q, r) = f.div_rem(&g)?;
                if r.is_zero() {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn check(&self, other: &PolyFp) -> PyResult<()> {
        if self.p != other.p {
            Err(PyValueError::new_err("Mixed base fields (different p)"))
        } else {
            Ok(())
        }
    }

    /// Return a compact debug-style representation.
    pub fn __repr__(&self) -> String {
        if self.coeffs.is_empty() {
            return format!("0 over F{}", self.p);
        }
        format!("{:?} over F{}", self.coeffs, self.p)
    }
}

fn add_coeffs(p: u64, a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut r = vec![0u64; n];
    for i in 0..n {
        let av = if i < a.len() { a[i] } else { 0 };
        let bv = if i < b.len() { b[i] } else { 0 };
        r[i] = (av + bv) % p;
    }
    r
}

fn mul_coeffs(p: u64, a: &[u64], b: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            let add = ((a[i] as u128 * b[j] as u128) % (p as u128)) as u64;
            let cur = r[i + j];
            r[i + j] = (cur + add) % p;
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_egcd_satisfies_bezout_identity() {
        let a = PolyFp::new(2, vec![1, 0, 1]).unwrap(); // 1 + x^2
        let b = PolyFp::new(2, vec![1, 1, 0, 1]).unwrap(); // 1 + x + x^3
        let (g, s, t) = a.egcd(&b).unwrap();

        let lhs = s
            .mul(&a)
            .unwrap()
            .add(&t.mul(&b).unwrap())
            .unwrap()
            .monic()
            .unwrap();
        assert_eq!(lhs, g.monic().unwrap());
        assert_eq!(g.monic().unwrap().coeffs(), vec![1]);
    }

    #[test]
    fn poly_gcd_matches_expected_common_factor() {
        let a = PolyFp::new(2, vec![0, 1, 1]).unwrap(); // x + x^2
        let b = PolyFp::new(2, vec![0, 1]).unwrap(); // x
        assert_eq!(a.gcd(&b).unwrap().coeffs(), vec![0, 1]);
    }

    #[test]
    fn poly_div_rem_identity_holds_over_small_inputs() {
        for p in [2_u64, 3] {
            for a0 in 0..p {
                for a1 in 0..p {
                    for a2 in 0..p {
                        for a3 in 0..p {
                            let a = PolyFp::new(
                                p,
                                vec![a0 as i128, a1 as i128, a2 as i128, a3 as i128],
                            )
                            .unwrap();
                            for b0 in 0..p {
                                for b1 in 0..p {
                                    for b2 in 0..p {
                                        let b = PolyFp::new(
                                            p,
                                            vec![b0 as i128, b1 as i128, b2 as i128],
                                        )
                                        .unwrap();
                                        if b.is_zero() {
                                            continue;
                                        }
                                        let (q, r) = a.div_rem(&b).unwrap();
                                        let recomposed = q.mul(&b).unwrap().add(&r).unwrap();
                                        assert_eq!(recomposed, a, "p={p}, a={:?}, b={:?}", a, b);
                                        assert!(r.is_zero() || r.degree() < b.degree());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

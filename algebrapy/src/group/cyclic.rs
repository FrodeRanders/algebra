use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn divisors_u64(n: u64) -> Vec<u64> {
    let mut small = Vec::new();
    let mut large = Vec::new();
    let mut d = 1u64;
    while d * d <= n {
        if n % d == 0 {
            small.push(d);
            if d * d != n {
                large.push(n / d);
            }
        }
        d += 1;
    }
    large.reverse();
    small.extend(large);
    small
}

/// The cyclic group `C_n` of order `n`, represented as `{0, ..., n-1}` under addition modulo `n`.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct Cn {
    n: u64,
}

#[pymethods]
impl Cn {
    #[new]
    pub fn new(n: u64) -> PyResult<Self> {
        if n < 1 {
            return Err(PyValueError::new_err("group order must be >= 1"));
        }
        Ok(Self { n })
    }

    pub fn order(&self) -> u64 {
        self.n
    }

    pub fn identity(&self) -> u64 {
        0
    }

    pub fn elements(&self) -> Vec<u64> {
        (0..self.n).collect()
    }

    pub fn generators(&self) -> Vec<u64> {
        (1..self.n).filter(|&k| gcd_u64(k, self.n) == 1).collect()
    }

    pub fn is_generator(&self, k: u64) -> bool {
        if k >= self.n {
            return false;
        }
        gcd_u64(k, self.n) == 1
    }

    pub fn element_order(&self, k: u64) -> PyResult<u64> {
        if k >= self.n {
            return Err(PyValueError::new_err("element out of range"));
        }
        let g = gcd_u64(k, self.n);
        Ok(self.n / g)
    }

    pub fn op(&self, a: u64, b: u64) -> PyResult<u64> {
        if a >= self.n || b >= self.n {
            return Err(PyValueError::new_err("element out of range"));
        }
        Ok((a + b) % self.n)
    }

    pub fn inv(&self, a: u64) -> PyResult<u64> {
        if a >= self.n {
            return Err(PyValueError::new_err("element out of range"));
        }
        Ok(if a == 0 { 0 } else { self.n - a })
    }

    pub fn pow(&self, a: u64, e: i64) -> PyResult<u64> {
        if a >= self.n {
            return Err(PyValueError::new_err("element out of range"));
        }
        let mut exp = e;
        if exp < 0 {
            exp = exp.rem_euclid(self.n as i64);
        }
        Ok(((a as u128 * exp as u128) % self.n as u128) as u64)
    }

    pub fn subgroup(&self, generators: Vec<u64>) -> PyResult<(Vec<u64>, Vec<u64>)> {
        if generators.is_empty() {
            return Ok((vec![0], vec![0]));
        }
        let d = generators
            .iter()
            .fold(0u64, |acc, &k| gcd_u64(acc, gcd_u64(k, self.n)));
        if d == 0 {
            return Ok((vec![0], vec![0]));
        }
        let elems: Vec<u64> = (0..self.n).step_by(d as usize).collect();
        Ok((vec![d % self.n], elems))
    }

    pub fn subgroups(&self) -> Vec<(u64, Vec<u64>)> {
        divisors_u64(self.n)
            .into_iter()
            .map(|d| {
                let g = if d == self.n { 0 } else { d };
                let elems: Vec<u64> = if d == 0 {
                    vec![0]
                } else {
                    (0..self.n).step_by(d as usize).collect()
                };
                (g, elems)
            })
            .collect()
    }

    pub fn __len__(&self) -> usize {
        self.n as usize
    }

    pub fn __repr__(&self) -> String {
        format!("C{}", self.n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cn12_generators_are_phi12() {
        let c = Cn::new(12).unwrap();
        let gens = c.generators();
        assert_eq!(gens, vec![1, 5, 7, 11]);
        assert!(c.is_generator(5));
        assert!(!c.is_generator(4));
    }

    #[test]
    fn cn12_element_order() {
        let c = Cn::new(12).unwrap();
        assert_eq!(c.element_order(0).unwrap(), 1);
        assert_eq!(c.element_order(4).unwrap(), 3);
        assert_eq!(c.element_order(6).unwrap(), 2);
        assert_eq!(c.element_order(5).unwrap(), 12);
    }

    #[test]
    fn cn12_subgroup_by_generator() {
        let c = Cn::new(12).unwrap();
        let (gens, elems) = c.subgroup(vec![8]).unwrap();
        assert_eq!(gens, vec![4]);
        assert_eq!(elems, vec![0, 4, 8]);
    }

    #[test]
    fn cn12_subgroups() {
        let c = Cn::new(12).unwrap();
        let subs = c.subgroups();
        assert_eq!(subs.len(), 6);
        let sizes: Vec<usize> = subs.iter().map(|(_g, e)| e.len()).collect();
        assert_eq!(sizes, vec![12, 6, 4, 3, 2, 1]);
    }

    #[test]
    fn cn1_trivial() {
        let c = Cn::new(1).unwrap();
        assert_eq!(c.elements(), vec![0]);
        assert_eq!(c.generators(), vec![] as Vec<u64>);
    }

    #[test]
    fn cn_ops() {
        let c = Cn::new(12).unwrap();
        assert_eq!(c.op(7, 8).unwrap(), 3);
        assert_eq!(c.inv(5).unwrap(), 7);
        assert_eq!(c.pow(5, 3).unwrap(), 3);
        assert_eq!(c.pow(5, -1).unwrap(), 7);
    }
}

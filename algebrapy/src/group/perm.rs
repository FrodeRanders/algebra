use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use std::collections::HashSet;

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Perm {
    n: usize,
    images: Vec<usize>,
}

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct Sn {
    n: usize,
}

fn is_bijection(n: usize, images: &[usize]) -> bool {
    if images.len() != n { return false; }
    let mut seen = vec![false; n];
    for &x in images {
        if x >= n || seen[x] { return false; }
        seen[x] = true;
    }
    true
}

#[pymethods]
impl Perm {
    #[new]
    pub fn new(n: usize, images: Vec<usize>) -> PyResult<Self> {
        if !is_bijection(n, &images) {
            return Err(PyValueError::new_err("images must be a bijection on 0..n-1"));
        }
        Ok(Self { n, images })
    }

    #[staticmethod]
    pub fn identity(n: usize) -> Self {
        Self { n, images: (0..n).collect() }
    }

    /// cycle given as list like [0,1,2] meaning (0 1 2)
    #[staticmethod]
    pub fn cycle(n: usize, cycle: Vec<usize>) -> PyResult<Self> {
        if cycle.len() < 2 {
            return Ok(Self::identity(n));
        }
        for &c in &cycle {
            if c >= n { return Err(PyValueError::new_err("cycle element out of range")); }
        }
        let mut images: Vec<usize> = (0..n).collect();
        for i in 0..cycle.len() {
            let a = cycle[i];
            let b = cycle[(i + 1) % cycle.len()];
            images[a] = b;
        }
        Ok(Self { n, images })
    }

    pub fn n(&self) -> usize { self.n }

    pub fn as_images(&self) -> Vec<usize> { self.images.clone() }

    /// Composition p ∘ q: apply q then p. (p.compose(q))(i) = p(q(i))
    pub fn compose(&self, other: &Perm) -> PyResult<Perm> {
        if self.n != other.n { return Err(PyValueError::new_err("different n")); }
        let n = self.n;
        let mut img = vec![0usize; n];
        for i in 0..n {
            img[i] = self.images[other.images[i]];
        }
        Ok(Perm { n, images: img })
    }

    pub fn inv(&self) -> Perm {
        let n = self.n;
        let mut inv = vec![0usize; n];
        for i in 0..n {
            inv[self.images[i]] = i;
        }
        Perm { n, images: inv }
    }

    pub fn pow(&self, exp: i64) -> PyResult<Perm> {
        if exp == 0 { return Ok(Perm::identity(self.n)); }
        if exp < 0 { return self.inv().pow(-exp); }
        let mut e = exp as u64;
        let mut result = Perm::identity(self.n);
        let mut base = self.clone();
        while e > 0 {
            if (e & 1) == 1 {
                result = result.compose(&base)?;
            }
            e >>= 1;
            if e > 0 {
                base = base.compose(&base)?;
            }
        }
        Ok(result)
    }

    pub fn order(&self) -> u64 {
        // lcm of cycle lengths
        let mut visited = vec![false; self.n];
        let mut l: u64 = 1;
        for i in 0..self.n {
            if visited[i] { continue; }
            let mut len = 0u64;
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                j = self.images[j];
                len += 1;
            }
            l = lcm(l, len);
        }
        l
    }

    pub fn __repr__(&self) -> String {
        format!("Perm(n={}, images={:?})", self.n, self.images)
    }

    // minimal dunders for convenience
    pub fn __mul__(&self, other: &Perm) -> PyResult<Perm> { self.compose(other) } // p*q = p∘q
}

#[pymethods]
impl Sn {
    #[new]
    pub fn new(n: usize) -> PyResult<Self> {
        if n > 8 {
            return Err(PyValueError::new_err("n too large for brute-force Sn enumeration (use n<=8)"));
        }
        Ok(Self { n })
    }

    pub fn n(&self) -> usize { self.n }

    pub fn identity(&self) -> Perm { Perm::identity(self.n) }

    pub fn elements(&self, max_size: Option<usize>) -> PyResult<Vec<Perm>> {
        let max = max_size.unwrap_or(5040);
        let mut a: Vec<usize> = (0..self.n).collect();
        let mut out = vec![Perm { n: self.n, images: a.clone() }];
        while next_permutation(&mut a) {
            if out.len() >= max {
                return Err(PyValueError::new_err("too many elements (increase max_size)"));
            }
            out.push(Perm { n: self.n, images: a.clone() });
        }
        Ok(out)
    }

    pub fn generated(&self, gens: Vec<Perm>) -> PyResult<Vec<Perm>> {
        self.generated_with_limit(gens, 5040)
    }

    pub fn generated_with_limit(&self, gens: Vec<Perm>, max_size: usize) -> PyResult<Vec<Perm>> {
        let max = max_size;
        for g in &gens {
            if g.n != self.n { return Err(PyValueError::new_err("generator has different n")); }
        }
        let id = self.identity();
        let mut set: HashSet<Perm> = HashSet::new();
        set.insert(id.clone());
        for g in &gens { set.insert(g.clone()); set.insert(g.inv()); }

        let mut changed = true;
        while changed {
            changed = false;
            let elems: Vec<Perm> = set.iter().cloned().collect();
            for a in &elems {
                for g in &gens {
                    let ag = a.compose(g)?;
                    if set.insert(ag) { changed = true; }
                    let ga = g.compose(a)?;
                    if set.insert(ga) { changed = true; }
                }
                if set.len() > max {
                    return Err(PyValueError::new_err("generated subgroup exceeds max_size"));
                }
            }
        }
        Ok(set.into_iter().collect())
    }

    pub fn __repr__(&self) -> String {
        format!("S{}", self.n)
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 { let t = a % b; a = b; b = t; }
    a
}
fn lcm(a: u64, b: u64) -> u64 {
    if a == 0 || b == 0 { 0 } else { (a / gcd(a,b)) * b }
}

// lexicographic next_permutation for Vec<usize>
fn next_permutation(a: &mut [usize]) -> bool {
    let n = a.len();
    if n < 2 { return false; }
    let mut i = n - 2;
    while a[i] >= a[i+1] {
        if i == 0 { a.reverse(); return false; }
        i -= 1;
    }
    let mut j = n - 1;
    while a[j] <= a[i] { j -= 1; }
    a.swap(i, j);
    a[i+1..].reverse();
    true
}

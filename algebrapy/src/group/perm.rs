use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use std::collections::HashSet;

/// A permutation of `{0, ..., n-1}` stored in image form.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Perm {
    n: usize,
    images: Vec<usize>,
}

/// The symmetric group `S_n` with small-scale enumeration and subgroup generation.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct Sn {
    n: usize,
}

fn is_bijection(n: usize, images: &[usize]) -> bool {
    if images.len() != n {
        return false;
    }
    let mut seen = vec![false; n];
    for &x in images {
        if x >= n || seen[x] {
            return false;
        }
        seen[x] = true;
    }
    true
}

#[pymethods]
impl Perm {
    #[new]
    /// Create a permutation from its image list.
    pub fn new(n: usize, images: Vec<usize>) -> PyResult<Self> {
        if !is_bijection(n, &images) {
            return Err(PyValueError::new_err(
                "images must be a bijection on 0..n-1",
            ));
        }
        Ok(Self { n, images })
    }

    #[staticmethod]
    /// Return the identity permutation on `n` points.
    pub fn identity(n: usize) -> Self {
        Self {
            n,
            images: (0..n).collect(),
        }
    }

    /// Build a cycle such as `[0, 1, 2]`, meaning `(0 1 2)`.
    #[staticmethod]
    pub fn cycle(n: usize, cycle: Vec<usize>) -> PyResult<Self> {
        if cycle.len() < 2 {
            return Ok(Self::identity(n));
        }
        for &c in &cycle {
            if c >= n {
                return Err(PyValueError::new_err("cycle element out of range"));
            }
        }
        let mut images: Vec<usize> = (0..n).collect();
        for i in 0..cycle.len() {
            let a = cycle[i];
            let b = cycle[(i + 1) % cycle.len()];
            images[a] = b;
        }
        Ok(Self { n, images })
    }

    /// Return the permutation degree.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Return the image list.
    pub fn as_images(&self) -> Vec<usize> {
        self.images.clone()
    }

    /// Composition p ∘ q: apply q then p. (p.compose(q))(i) = p(q(i))
    pub fn compose(&self, other: &Perm) -> PyResult<Perm> {
        if self.n != other.n {
            return Err(PyValueError::new_err("different n"));
        }
        let n = self.n;
        let mut img = vec![0usize; n];
        for i in 0..n {
            img[i] = self.images[other.images[i]];
        }
        Ok(Perm { n, images: img })
    }

    /// Return the inverse permutation.
    pub fn inv(&self) -> Perm {
        let n = self.n;
        let mut inv = vec![0usize; n];
        for i in 0..n {
            inv[self.images[i]] = i;
        }
        Perm { n, images: inv }
    }

    /// Raise the permutation to an integer power.
    pub fn pow(&self, exp: i64) -> PyResult<Perm> {
        if exp == 0 {
            return Ok(Perm::identity(self.n));
        }
        if exp < 0 {
            return self.inv().pow(-exp);
        }
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

    /// Return the conjugate `g * self * g^-1`.
    pub fn conjugate_by(&self, g: &Perm) -> PyResult<Perm> {
        if self.n != g.n {
            return Err(PyValueError::new_err("different n"));
        }
        g.compose(self)?.compose(&g.inv())
    }

    /// Return the order of the permutation.
    pub fn order(&self) -> u64 {
        // lcm of cycle lengths
        let mut visited = vec![false; self.n];
        let mut l: u64 = 1;
        for i in 0..self.n {
            if visited[i] {
                continue;
            }
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

    /// Return the nontrivial disjoint cycles of the permutation.
    pub fn cycles(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.n];
        let mut out: Vec<Vec<usize>> = Vec::new();
        for i in 0..self.n {
            if visited[i] {
                continue;
            }
            let mut cycle: Vec<usize> = Vec::new();
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                cycle.push(j);
                j = self.images[j];
            }
            if cycle.len() > 1 {
                out.push(cycle);
            }
        }
        out
    }

    /// Return the sorted lengths of the nontrivial disjoint cycles.
    pub fn cycle_type(&self) -> Vec<usize> {
        let mut ty: Vec<usize> = self.cycles().into_iter().map(|c| c.len()).collect();
        ty.sort_unstable();
        ty
    }

    /// Return the permutation in cycle notation, omitting fixed points.
    pub fn cycle_notation(&self) -> String {
        let cycles = self.cycles();
        if cycles.is_empty() {
            return "()".to_string();
        }
        cycles
            .into_iter()
            .map(|cycle| {
                format!(
                    "({})",
                    cycle
                        .into_iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(" ")
                )
            })
            .collect::<Vec<_>>()
            .join("")
    }

    /// Return a debug-style representation.
    pub fn __repr__(&self) -> String {
        format!("Perm(n={}, images={:?})", self.n, self.images)
    }

    /// Compose two permutations with `p * q = p.compose(q)`.
    pub fn __mul__(&self, other: &Perm) -> PyResult<Perm> {
        self.compose(other)
    } // p*q = p∘q
}

#[pymethods]
impl Sn {
    #[new]
    /// Construct the symmetric group `S_n`.
    pub fn new(n: usize) -> PyResult<Self> {
        Ok(Self { n })
    }

    /// Return the degree `n`.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Return the identity element of `S_n`.
    pub fn identity(&self) -> Perm {
        Perm::identity(self.n)
    }

    /// Enumerate all elements of `S_n`, subject to `max_size`.
    pub fn elements(&self, max_size: Option<usize>) -> PyResult<Vec<Perm>> {
        let max = max_size.unwrap_or(5040);
        let mut a: Vec<usize> = (0..self.n).collect();
        let mut out = vec![Perm {
            n: self.n,
            images: a.clone(),
        }];
        while next_permutation(&mut a) {
            if out.len() >= max {
                return Err(PyValueError::new_err(
                    "too many elements (increase max_size)",
                ));
            }
            out.push(Perm {
                n: self.n,
                images: a.clone(),
            });
        }
        Ok(out)
    }

    /// Return the subgroup generated by `gens` with a default size bound.
    pub fn generated(&self, gens: Vec<Perm>) -> PyResult<Vec<Perm>> {
        self.generated_with_limit(gens, 5040)
    }

    /// Return the subgroup generated by `gens`, failing if it exceeds `max_size`.
    pub fn generated_with_limit(&self, gens: Vec<Perm>, max_size: usize) -> PyResult<Vec<Perm>> {
        let max = max_size;
        for g in &gens {
            if g.n != self.n {
                return Err(PyValueError::new_err("generator has different n"));
            }
        }
        let id = self.identity();
        let mut set: HashSet<Perm> = HashSet::new();
        let mut queue: Vec<Perm> = Vec::new();
        let mut closure_gens: Vec<Perm> = Vec::new();

        set.insert(id.clone());
        queue.push(id);
        for g in &gens {
            let g_inv = g.inv();
            closure_gens.push(g.clone());
            closure_gens.push(g_inv.clone());
            if set.insert(g.clone()) {
                queue.push(g.clone());
            }
            if set.insert(g_inv.clone()) {
                queue.push(g_inv);
            }
        }

        while let Some(a) = queue.pop() {
            for g in &closure_gens {
                let ag = a.compose(g)?;
                if set.insert(ag.clone()) {
                    if set.len() > max {
                        return Err(PyValueError::new_err("generated subgroup exceeds max_size"));
                    }
                    queue.push(ag);
                }
                let ga = g.compose(&a)?;
                if set.insert(ga.clone()) {
                    if set.len() > max {
                        return Err(PyValueError::new_err("generated subgroup exceeds max_size"));
                    }
                    queue.push(ga);
                }
            }
        }
        Ok(set.into_iter().collect())
    }

    /// Return the orbit of `point` under the subgroup generated by `gens`.
    pub fn orbit(&self, point: usize, gens: Vec<Perm>) -> PyResult<Vec<usize>> {
        if point >= self.n {
            return Err(PyValueError::new_err("point out of range"));
        }
        let subgroup = self.generated(gens)?;
        let mut orbit: Vec<usize> = subgroup.into_iter().map(|g| g.images[point]).collect();
        orbit.sort_unstable();
        orbit.dedup();
        Ok(orbit)
    }

    /// Return the orbits of the subgroup generated by `gens` on `{0, ..., n-1}`.
    pub fn orbits(&self, gens: Vec<Perm>) -> PyResult<Vec<Vec<usize>>> {
        let subgroup = self.generated(gens)?;
        let mut visited = vec![false; self.n];
        let mut out: Vec<Vec<usize>> = Vec::new();

        for start in 0..self.n {
            if visited[start] {
                continue;
            }
            let mut orbit: Vec<usize> = subgroup.iter().map(|g| g.images[start]).collect();
            orbit.sort_unstable();
            orbit.dedup();
            for &x in &orbit {
                visited[x] = true;
            }
            out.push(orbit);
        }

        Ok(out)
    }

    /// Return the stabilizer subgroup of `point` inside the subgroup generated by `gens`.
    pub fn stabilizer(&self, point: usize, gens: Vec<Perm>) -> PyResult<Vec<Perm>> {
        if point >= self.n {
            return Err(PyValueError::new_err("point out of range"));
        }
        let subgroup = self.generated(gens)?;
        Ok(subgroup
            .into_iter()
            .filter(|g| g.images[point] == point)
            .collect())
    }

    /// Return the size of the stabilizer of `point` in the subgroup generated by `gens`.
    pub fn stabilizer_size(&self, point: usize, gens: Vec<Perm>) -> PyResult<usize> {
        Ok(self.stabilizer(point, gens)?.len())
    }

    /// Return the conjugacy class of `perm` inside `S_n`, subject to `max_size`.
    pub fn conjugacy_class(&self, perm: &Perm, max_size: Option<usize>) -> PyResult<Vec<Perm>> {
        if perm.n != self.n {
            return Err(PyValueError::new_err("permutation has different n"));
        }
        let mut class: HashSet<Perm> = HashSet::new();
        for g in self.elements(max_size)? {
            class.insert(perm.conjugate_by(&g)?);
        }
        Ok(class.into_iter().collect())
    }

    /// Return the size of the conjugacy class of `perm` inside `S_n`, subject to `max_size`.
    pub fn conjugacy_class_size(&self, perm: &Perm, max_size: Option<usize>) -> PyResult<usize> {
        Ok(self.conjugacy_class(perm, max_size)?.len())
    }

    /// Return the size of the subgroup generated by `gens`.
    pub fn subgroup_size(&self, gens: Vec<Perm>) -> PyResult<usize> {
        Ok(self.generated(gens)?.len())
    }

    /// Return whether the subgroup generated by `gens` acts transitively.
    pub fn is_transitive(&self, gens: Vec<Perm>) -> PyResult<bool> {
        Ok(self.orbits(gens)?.len() == 1)
    }

    /// Return a short representation such as `S4`.
    pub fn __repr__(&self) -> String {
        format!("S{}", self.n)
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}
fn lcm(a: u64, b: u64) -> u64 {
    if a == 0 || b == 0 {
        0
    } else {
        (a / gcd(a, b)) * b
    }
}

// lexicographic next_permutation for Vec<usize>
fn next_permutation(a: &mut [usize]) -> bool {
    let n = a.len();
    if n < 2 {
        return false;
    }
    let mut i = n - 2;
    while a[i] >= a[i + 1] {
        if i == 0 {
            a.reverse();
            return false;
        }
        i -= 1;
    }
    let mut j = n - 1;
    while a[j] <= a[i] {
        j -= 1;
    }
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generated_subgroup_contains_inverse_words() {
        let s4 = Sn::new(4).unwrap();
        let a = Perm::cycle(4, vec![0, 1, 2]).unwrap();
        let b = Perm::cycle(4, vec![1, 2, 3]).unwrap();

        let subgroup = s4.generated(vec![a.clone(), b.clone()]).unwrap();
        let subgroup_set: HashSet<Perm> = subgroup.into_iter().collect();

        let target = a.inv().compose(&b).unwrap().compose(&a.inv()).unwrap();
        assert!(subgroup_set.contains(&target));
    }

    #[test]
    fn generated_subgroup_of_adjacent_transpositions_is_s4() {
        let s4 = Sn::new(4).unwrap();
        let t01 = Perm::cycle(4, vec![0, 1]).unwrap();
        let t12 = Perm::cycle(4, vec![1, 2]).unwrap();
        let t23 = Perm::cycle(4, vec![2, 3]).unwrap();

        let subgroup = s4.generated(vec![t01, t12, t23]).unwrap();
        assert_eq!(subgroup.len(), 24);
    }

    #[test]
    fn orbits_of_generated_action_are_computed() {
        let s6 = Sn::new(6).unwrap();
        let a = Perm::cycle(6, vec![0, 1, 2]).unwrap();
        let b = Perm::cycle(6, vec![3, 4]).unwrap();

        let mut orbits = s6.orbits(vec![a, b]).unwrap();
        orbits.sort();
        assert_eq!(orbits, vec![vec![0, 1, 2], vec![3, 4], vec![5]]);
    }

    #[test]
    fn orbit_of_point_is_computed() {
        let s5 = Sn::new(5).unwrap();
        let g = Perm::cycle(5, vec![1, 3, 4]).unwrap();
        assert_eq!(s5.orbit(1, vec![g]).unwrap(), vec![1, 3, 4]);
    }

    #[test]
    fn subgroup_size_and_transitivity_are_reported() {
        let s4 = Sn::new(4).unwrap();
        let g = Perm::cycle(4, vec![0, 1, 2, 3]).unwrap();
        assert_eq!(s4.subgroup_size(vec![g.clone()]).unwrap(), 4);
        assert!(s4.is_transitive(vec![g]).unwrap());
    }

    #[test]
    fn stabilizer_is_computed_and_matches_orbit_stabilizer() {
        let s4 = Sn::new(4).unwrap();
        let a = Perm::cycle(4, vec![0, 1, 2, 3]).unwrap();
        let b = Perm::cycle(4, vec![1, 3]).unwrap();
        let gens = vec![a.clone(), b.clone()];

        let subgroup_size = s4.subgroup_size(gens.clone()).unwrap();
        let orbit = s4.orbit(0, gens.clone()).unwrap();
        let stabilizer = s4.stabilizer(0, gens.clone()).unwrap();

        assert_eq!(stabilizer.len(), s4.stabilizer_size(0, gens).unwrap());
        assert_eq!(subgroup_size, orbit.len() * stabilizer.len());
    }

    #[test]
    fn cycle_structure_ignores_fixed_points() {
        let p = Perm::new(6, vec![1, 0, 2, 4, 5, 3]).unwrap();
        assert_eq!(p.cycles(), vec![vec![0, 1], vec![3, 4, 5]]);
        assert_eq!(p.cycle_type(), vec![2, 3]);
        assert_eq!(p.order(), 6);
    }

    #[test]
    fn cycle_notation_omits_fixed_points() {
        let p = Perm::new(6, vec![1, 0, 2, 4, 5, 3]).unwrap();
        assert_eq!(p.cycle_notation(), "(0 1)(3 4 5)");
        assert_eq!(Perm::identity(4).cycle_notation(), "()");
    }

    #[test]
    fn conjugation_preserves_cycle_type() {
        let p = Perm::cycle(5, vec![0, 1, 2]).unwrap();
        let g = Perm::cycle(5, vec![0, 4, 3, 2, 1]).unwrap();
        assert_eq!(p.conjugate_by(&g).unwrap().cycle_type(), p.cycle_type());
    }

    #[test]
    fn conjugacy_class_size_of_transposition_in_s4_is_six() {
        let s4 = Sn::new(4).unwrap();
        let t = Perm::cycle(4, vec![0, 1]).unwrap();
        assert_eq!(s4.conjugacy_class_size(&t, Some(24)).unwrap(), 6);
    }
}

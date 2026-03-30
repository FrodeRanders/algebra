use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use std::collections::HashSet;

use crate::arith::prime::is_prime_u64;

/// A permutation of `{0, ..., n-1}` stored in image form.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Perm {
    n: usize,
    images: Vec<usize>,
}

/// An explicit finite subgroup of `S_n`.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct PermSubgroup {
    n: usize,
    elements: Vec<Perm>,
    generators: Vec<Perm>,
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

fn perm_key(perm: &Perm) -> Vec<usize> {
    perm.images.clone()
}

fn subgroup_key(elements: &[Perm]) -> Vec<Vec<usize>> {
    let mut key: Vec<Vec<usize>> = elements.iter().map(perm_key).collect();
    key.sort_unstable();
    key
}

fn generated_elements(n: usize, gens: &[Perm], max_size: usize) -> PyResult<Vec<Perm>> {
    for g in gens {
        if g.n != n {
            return Err(PyValueError::new_err("generator has different n"));
        }
    }

    let id = Perm::identity(n);
    let mut set: HashSet<Perm> = HashSet::new();
    let mut queue: Vec<Perm> = Vec::new();
    let mut closure_gens: Vec<Perm> = Vec::new();

    set.insert(id.clone());
    queue.push(id);
    for g in gens {
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
                if set.len() > max_size {
                    return Err(PyValueError::new_err("generated subgroup exceeds max_size"));
                }
                queue.push(ag);
            }
            let ga = g.compose(&a)?;
            if set.insert(ga.clone()) {
                if set.len() > max_size {
                    return Err(PyValueError::new_err("generated subgroup exceeds max_size"));
                }
                queue.push(ga);
            }
        }
    }

    let mut elements: Vec<Perm> = set.into_iter().collect();
    elements.sort_unstable_by_key(perm_key);
    Ok(elements)
}

fn validate_p(p: u64) -> PyResult<()> {
    if !is_prime_u64(p) {
        return Err(PyValueError::new_err("p must be prime"));
    }
    Ok(())
}

fn p_part(mut n: usize, p: u64) -> usize {
    let p_usize = p as usize;
    let mut out = 1usize;
    while n % p_usize == 0 {
        out *= p_usize;
        n /= p_usize;
    }
    out
}

impl PermSubgroup {
    fn from_generated(n: usize, gens: Vec<Perm>, max_size: usize) -> PyResult<Self> {
        let elements = generated_elements(n, &gens, max_size)?;
        Ok(Self {
            n,
            elements,
            generators: gens,
        })
    }

    fn from_elements(n: usize, elements: Vec<Perm>, generators: Vec<Perm>) -> PyResult<Self> {
        if elements.is_empty() {
            return Err(PyValueError::new_err("subgroup must contain at least one element"));
        }
        for g in &elements {
            if g.n != n {
                return Err(PyValueError::new_err("element has different n"));
            }
        }
        let mut sorted = elements;
        sorted.sort_unstable_by_key(perm_key);
        sorted.dedup();
        let regen = generated_elements(n, &sorted, sorted.len())?;
        if subgroup_key(&regen) != subgroup_key(&sorted) {
            return Err(PyValueError::new_err("elements do not form a subgroup"));
        }
        Ok(Self {
            n,
            elements: sorted,
            generators,
        })
    }

    fn contains_perm(&self, perm: &Perm) -> bool {
        self.elements.contains(perm)
    }

    fn ensure_same_degree(&self, other: &PermSubgroup) -> PyResult<()> {
        if self.n != other.n {
            return Err(PyValueError::new_err("different n"));
        }
        Ok(())
    }

    fn ensure_subset_of(&self, parent: &PermSubgroup) -> PyResult<()> {
        self.ensure_same_degree(parent)?;
        if self.elements.iter().all(|g| parent.contains_perm(g)) {
            Ok(())
        } else {
            Err(PyValueError::new_err("subgroup is not contained in parent"))
        }
    }

    fn trivial(n: usize) -> Self {
        let id = Perm::identity(n);
        Self {
            n,
            elements: vec![id.clone()],
            generators: vec![id],
        }
    }

    fn enumerate_subgroups(&self, max_subgroups: usize) -> PyResult<Vec<PermSubgroup>> {
        let mut seen: HashSet<Vec<Vec<usize>>> = HashSet::new();
        let mut out: Vec<PermSubgroup> = Vec::new();

        let trivial = PermSubgroup::trivial(self.n);
        seen.insert(subgroup_key(&trivial.elements));
        out.push(trivial);

        let mut idx = 0usize;
        while idx < out.len() {
            let current = out[idx].clone();
            for g in &self.elements {
                if current.contains_perm(g) {
                    continue;
                }
                let mut gens = current.elements.clone();
                gens.push(g.clone());
                let candidate = PermSubgroup::from_generated(self.n, gens.clone(), self.order())?;
                let key = subgroup_key(&candidate.elements);
                if seen.insert(key) {
                    out.push(PermSubgroup {
                        n: self.n,
                        elements: candidate.elements,
                        generators: gens,
                    });
                    if out.len() > max_subgroups {
                        return Err(PyValueError::new_err(
                            "subgroup enumeration exceeds max_size",
                        ));
                    }
                }
            }
            idx += 1;
        }

        Ok(out)
    }
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
        for (i, slot) in img.iter_mut().enumerate().take(n) {
            *slot = self.images[other.images[i]];
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
    }
}

#[pymethods]
impl PermSubgroup {
    /// Return the ambient permutation degree.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Return the explicit subgroup elements.
    pub fn elements(&self) -> Vec<Perm> {
        self.elements.clone()
    }

    /// Return the generators used to construct the subgroup when known.
    pub fn generators(&self) -> Vec<Perm> {
        self.generators.clone()
    }

    /// Return the subgroup order.
    pub fn order(&self) -> usize {
        self.elements.len()
    }

    /// Return whether `perm` lies in the subgroup.
    pub fn contains(&self, perm: &Perm) -> PyResult<bool> {
        if perm.n != self.n {
            return Err(PyValueError::new_err("different n"));
        }
        Ok(self.contains_perm(perm))
    }

    /// Return whether the subgroup is abelian.
    pub fn is_abelian(&self) -> PyResult<bool> {
        for a in &self.elements {
            for b in &self.elements {
                if a.compose(b)? != b.compose(a)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Return the orbit of `point` under the subgroup action.
    pub fn orbit(&self, point: usize) -> PyResult<Vec<usize>> {
        if point >= self.n {
            return Err(PyValueError::new_err("point out of range"));
        }
        let mut orbit: Vec<usize> = self.elements.iter().map(|g| g.images[point]).collect();
        orbit.sort_unstable();
        orbit.dedup();
        Ok(orbit)
    }

    /// Return the orbits of the subgroup action on `{0, ..., n-1}`.
    pub fn orbits(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.n];
        let mut out: Vec<Vec<usize>> = Vec::new();

        for start in 0..self.n {
            if visited[start] {
                continue;
            }
            let mut orbit: Vec<usize> = self.elements.iter().map(|g| g.images[start]).collect();
            orbit.sort_unstable();
            orbit.dedup();
            for &x in &orbit {
                visited[x] = true;
            }
            out.push(orbit);
        }

        out
    }

    /// Return the stabilizer subgroup of `point`.
    pub fn stabilizer(&self, point: usize) -> PyResult<PermSubgroup> {
        if point >= self.n {
            return Err(PyValueError::new_err("point out of range"));
        }
        let elements: Vec<Perm> = self
            .elements
            .iter()
            .filter(|g| g.images[point] == point)
            .cloned()
            .collect();
        PermSubgroup::from_elements(self.n, elements.clone(), elements)
    }

    /// Return the size of the stabilizer of `point`.
    pub fn stabilizer_size(&self, point: usize) -> PyResult<usize> {
        Ok(self.stabilizer(point)?.order())
    }

    /// Return whether the subgroup action is transitive.
    pub fn is_transitive(&self) -> bool {
        self.orbits().len() == 1
    }

    /// Return the conjugate subgroup `g * H * g^-1`.
    pub fn conjugate_by(&self, g: &Perm) -> PyResult<PermSubgroup> {
        if g.n != self.n {
            return Err(PyValueError::new_err("different n"));
        }
        let elements: Vec<Perm> = self
            .elements
            .iter()
            .map(|h| h.conjugate_by(g))
            .collect::<PyResult<Vec<_>>>()?;
        let generators: Vec<Perm> = self
            .generators
            .iter()
            .map(|h| h.conjugate_by(g))
            .collect::<PyResult<Vec<_>>>()?;
        PermSubgroup::from_elements(self.n, elements, generators)
    }

    /// Return the intersection with `other`.
    pub fn intersection(&self, other: &PermSubgroup) -> PyResult<PermSubgroup> {
        self.ensure_same_degree(other)?;
        let elements: Vec<Perm> = self
            .elements
            .iter()
            .filter(|g| other.contains_perm(g))
            .cloned()
            .collect();
        PermSubgroup::from_elements(self.n, elements.clone(), elements)
    }

    /// Return whether the subgroup is normal in `parent`.
    pub fn is_normal_in(&self, parent: &PermSubgroup) -> PyResult<bool> {
        self.ensure_subset_of(parent)?;
        for g in &parent.elements {
            for h in &self.elements {
                let conj = h.conjugate_by(g)?;
                if !self.contains_perm(&conj) {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Return the highest power of `p` dividing the subgroup order.
    pub fn p_part_order(&self, p: u64) -> PyResult<usize> {
        validate_p(p)?;
        Ok(p_part(self.order(), p))
    }

    /// Return whether the subgroup order is a power of `p`.
    pub fn is_p_group(&self, p: u64) -> PyResult<bool> {
        validate_p(p)?;
        Ok(self.p_part_order(p)? == self.order())
    }

    /// Return whether `self` is a Sylow `p`-subgroup of `parent`.
    pub fn is_sylow_p_subgroup(&self, p: u64, parent: &PermSubgroup) -> PyResult<bool> {
        validate_p(p)?;
        self.ensure_subset_of(parent)?;
        Ok(self.is_p_group(p)? && self.order() == parent.p_part_order(p)?)
    }

    /// Enumerate the Sylow `p`-subgroups of this explicit finite subgroup.
    #[pyo3(signature = (p, max_size=None))]
    pub fn sylow_p_subgroups(
        &self,
        p: u64,
        max_size: Option<usize>,
    ) -> PyResult<Vec<PermSubgroup>> {
        validate_p(p)?;
        let target_order = self.p_part_order(p)?;
        let max_subgroups = max_size.unwrap_or(4096);
        let subgroups = self.enumerate_subgroups(max_subgroups)?;
        Ok(subgroups
            .into_iter()
            .filter(|h| h.order() == target_order)
            .collect())
    }

    pub fn __len__(&self) -> usize {
        self.order()
    }

    pub fn __repr__(&self) -> String {
        format!("PermSubgroup(n={}, order={})", self.n, self.order())
    }
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
    pub fn generated(&self, gens: Vec<Perm>) -> PyResult<PermSubgroup> {
        self.generated_with_limit(gens, 5040)
    }

    /// Return the subgroup generated by `gens`, failing if it exceeds `max_size`.
    pub fn generated_with_limit(&self, gens: Vec<Perm>, max_size: usize) -> PyResult<PermSubgroup> {
        PermSubgroup::from_generated(self.n, gens, max_size)
    }

    /// Return the orbit of `point` under the subgroup generated by `gens`.
    pub fn orbit(&self, point: usize, gens: Vec<Perm>) -> PyResult<Vec<usize>> {
        self.generated(gens)?.orbit(point)
    }

    /// Return the orbits of the subgroup generated by `gens` on `{0, ..., n-1}`.
    pub fn orbits(&self, gens: Vec<Perm>) -> PyResult<Vec<Vec<usize>>> {
        Ok(self.generated(gens)?.orbits())
    }

    /// Return the stabilizer subgroup of `point` inside the subgroup generated by `gens`.
    pub fn stabilizer(&self, point: usize, gens: Vec<Perm>) -> PyResult<PermSubgroup> {
        self.generated(gens)?.stabilizer(point)
    }

    /// Return the size of the stabilizer of `point` in the subgroup generated by `gens`.
    pub fn stabilizer_size(&self, point: usize, gens: Vec<Perm>) -> PyResult<usize> {
        self.generated(gens)?.stabilizer_size(point)
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
        Ok(self.generated(gens)?.order())
    }

    /// Return whether the subgroup generated by `gens` acts transitively.
    pub fn is_transitive(&self, gens: Vec<Perm>) -> PyResult<bool> {
        Ok(self.generated(gens)?.is_transitive())
    }

    /// Enumerate the Sylow `p`-subgroups of the subgroup generated by `gens`.
    #[pyo3(signature = (gens, p, max_size=None))]
    pub fn sylow_p_subgroups_of_generated(
        &self,
        gens: Vec<Perm>,
        p: u64,
        max_size: Option<usize>,
    ) -> PyResult<Vec<PermSubgroup>> {
        self.generated(gens)?.sylow_p_subgroups(p, max_size)
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

    fn dihedral_8_in_s4() -> PermSubgroup {
        let s4 = Sn::new(4).unwrap();
        let r = Perm::cycle(4, vec![0, 1, 2, 3]).unwrap();
        let f = Perm::new(4, vec![0, 3, 2, 1]).unwrap();
        s4.generated(vec![r, f]).unwrap()
    }

    #[test]
    fn generated_subgroup_contains_inverse_words() {
        let s4 = Sn::new(4).unwrap();
        let a = Perm::cycle(4, vec![0, 1, 2]).unwrap();
        let b = Perm::cycle(4, vec![1, 2, 3]).unwrap();

        let subgroup = s4.generated(vec![a.clone(), b.clone()]).unwrap();

        let target = a.inv().compose(&b).unwrap().compose(&a.inv()).unwrap();
        assert!(subgroup.contains(&target).unwrap());
    }

    #[test]
    fn generated_subgroup_of_adjacent_transpositions_is_s4() {
        let s4 = Sn::new(4).unwrap();
        let t01 = Perm::cycle(4, vec![0, 1]).unwrap();
        let t12 = Perm::cycle(4, vec![1, 2]).unwrap();
        let t23 = Perm::cycle(4, vec![2, 3]).unwrap();

        let subgroup = s4.generated(vec![t01, t12, t23]).unwrap();
        assert_eq!(subgroup.order(), 24);
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

        let group = s4.generated(gens.clone()).unwrap();
        let subgroup_size = group.order();
        let orbit = s4.orbit(0, gens.clone()).unwrap();
        let stabilizer = s4.stabilizer(0, gens).unwrap();

        assert_eq!(stabilizer.order(), s4.stabilizer_size(0, vec![a, b]).unwrap());
        assert_eq!(subgroup_size, orbit.len() * stabilizer.order());
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

    #[test]
    fn subgroup_reports_basic_structure() {
        let group = dihedral_8_in_s4();
        assert!(group.contains(&Perm::identity(4)).unwrap());
        assert!(!group.is_abelian().unwrap());
        assert_eq!(group.p_part_order(2).unwrap(), 8);
        assert!(group.is_p_group(2).unwrap());
        assert!(!group.is_p_group(3).unwrap());
    }

    #[test]
    fn subgroup_conjugation_and_intersection_work() {
        let s4 = Sn::new(4).unwrap();
        let h = s4.generated(vec![Perm::cycle(4, vec![0, 1]).unwrap()]).unwrap();
        let g = Perm::cycle(4, vec![0, 2]).unwrap();
        let conj = h.conjugate_by(&g).unwrap();
        let expected = Perm::cycle(4, vec![1, 2]).unwrap();

        assert!(conj.contains(&expected).unwrap());
        assert_eq!(h.intersection(&conj).unwrap().order(), 1);
    }

    #[test]
    fn normality_is_detected() {
        let s4 = Sn::new(4).unwrap();
        let v4 = s4
            .generated(vec![
                Perm::new(4, vec![1, 0, 3, 2]).unwrap(),
                Perm::new(4, vec![2, 3, 0, 1]).unwrap(),
            ])
            .unwrap();
        let s4_full = s4
            .generated(vec![
                Perm::cycle(4, vec![0, 1]).unwrap(),
                Perm::cycle(4, vec![0, 1, 2, 3]).unwrap(),
            ])
            .unwrap();

        assert!(v4.is_normal_in(&s4_full).unwrap());
    }

    #[test]
    fn sylow_subgroups_of_s3_are_found() {
        let s3 = Sn::new(3).unwrap();
        let sigma = Perm::cycle(3, vec![0, 1, 2]).unwrap();
        let tau = Perm::cycle(3, vec![0, 1]).unwrap();
        let g = s3.generated(vec![sigma, tau]).unwrap();

        let sylow2 = g.sylow_p_subgroups(2, None).unwrap();
        let sylow3 = g.sylow_p_subgroups(3, None).unwrap();

        assert_eq!(sylow2.len(), 3);
        assert!(sylow2.iter().all(|h| h.order() == 2));
        assert_eq!(sylow3.len(), 1);
        assert_eq!(sylow3[0].order(), 3);
        assert!(sylow3[0].is_normal_in(&g).unwrap());
    }

    #[test]
    fn sylow_subgroups_of_s4_have_expected_counts() {
        let s4 = Sn::new(4).unwrap();
        let g = s4
            .generated(vec![
                Perm::cycle(4, vec![0, 1]).unwrap(),
                Perm::cycle(4, vec![0, 1, 2, 3]).unwrap(),
            ])
            .unwrap();

        let sylow2 = g.sylow_p_subgroups(2, None).unwrap();
        let sylow3 = g.sylow_p_subgroups(3, None).unwrap();

        assert_eq!(sylow2.len(), 3);
        assert!(sylow2.iter().all(|h| h.order() == 8));
        assert_eq!(sylow3.len(), 4);
        assert!(sylow3.iter().all(|h| h.order() == 3));
        assert!(sylow2.iter().all(|h| h.is_sylow_p_subgroup(2, &g).unwrap()));
    }
}

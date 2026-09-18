use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::arith::prime::is_prime_u64;
use crate::group::perm::{Perm, PermSubgroup};

type Images = Vec<usize>;

#[derive(Clone, Debug)]
struct StabilizerLevel {
    point: usize,
    transversal: BTreeMap<usize, Images>,
}

fn identity_images(n: usize) -> Images {
    (0..n).collect()
}

fn inverse_images(g: &[usize]) -> Images {
    let mut inv = vec![0usize; g.len()];
    for (i, &x) in g.iter().enumerate() {
        inv[x] = i;
    }
    inv
}

/// `a ∘ b`, so `result[i] = a[b[i]]` (apply `b`, then `a`).
fn compose_images(a: &[usize], b: &[usize]) -> Images {
    b.iter().map(|&x| a[x]).collect()
}

fn fixes_pointwise(g: &[usize], points: &[usize]) -> bool {
    points.iter().all(|&p| g[p] == p)
}

fn first_moved_point(g: &[usize]) -> Option<usize> {
    (0..g.len()).find(|&i| g[i] != i)
}

fn orbit_transversal(
    n: usize,
    point: usize,
    gens: &[Images],
) -> (Vec<usize>, BTreeMap<usize, Images>) {
    let mut transversal: BTreeMap<usize, Images> = BTreeMap::new();
    transversal.insert(point, identity_images(n));
    let mut queue: VecDeque<usize> = VecDeque::new();
    queue.push_back(point);

    while let Some(x) = queue.pop_front() {
        let ux = transversal
            .get(&x)
            .expect("orbit point missing from transversal")
            .clone();
        for s in gens {
            let y = s[x];
            if let Entry::Vacant(slot) = transversal.entry(y) {
                slot.insert(compose_images(s, &ux));
                queue.push_back(y);
            }
        }
    }

    let orbit: Vec<usize> = transversal.keys().copied().collect();
    (orbit, transversal)
}

type ChainData = (Vec<usize>, Vec<Images>, Vec<StabilizerLevel>, u128);

fn build_perm_group(n: usize, generators: &[Perm]) -> Result<ChainData, String> {
    let id = identity_images(n);
    let mut strong: BTreeSet<Images> = BTreeSet::new();

    for g in generators {
        let images = g.as_images();
        if images != id {
            strong.insert(images.clone());
            strong.insert(inverse_images(&images));
        }
    }

    loop {
        // Build the stabilizer chain for the current generating set.
        let mut base: Vec<usize> = Vec::new();
        let mut levels: Vec<StabilizerLevel> = Vec::new();

        loop {
            let mut next_point: Option<usize> = None;
            for s in &strong {
                if fixes_pointwise(s, &base)
                    && let Some(point) = first_moved_point(s)
                {
                    next_point = Some(point);
                    break;
                }
            }
            let Some(point) = next_point else {
                break;
            };
            let level_gens: Vec<Images> = strong
                .iter()
                .filter(|s| fixes_pointwise(s, &base))
                .cloned()
                .collect();
            let (_orbit, transversal) = orbit_transversal(n, point, &level_gens);
            base.push(point);
            levels.push(StabilizerLevel { point, transversal });
        }

        // Sift every Schreier generator through the chain below its level.
        // The first generator that does not sift to the identity is a new
        // strong generator, so add it and rebuild.
        let mut new_generator: Option<Images> = None;
        let mut failure_level = (0usize, 0usize);

        'verify: for i in 0..base.len() {
            let level_gens: Vec<Images> = strong
                .iter()
                .filter(|s| fixes_pointwise(s, &base[..i]))
                .cloned()
                .collect();
            let orbit: Vec<usize> = levels[i].transversal.keys().copied().collect();

            for x in orbit {
                let ux = &levels[i].transversal[&x];
                for s in &level_gens {
                    let y = s[x];
                    let uy = &levels[i].transversal[&y];
                    let mut h = compose_images(&compose_images(&inverse_images(uy), s), ux);

                    let mut j = i + 1;
                    while j < base.len() {
                        let image = h[base[j]];
                        match levels[j].transversal.get(&image) {
                            Some(u) => {
                                h = compose_images(&inverse_images(u), &h);
                                j += 1;
                            }
                            None => break,
                        }
                    }

                    if j < base.len() || h != id {
                        new_generator = Some(h);
                        failure_level = (i, j);
                        break 'verify;
                    }
                }
            }
        }

        match new_generator {
            Some(h) => {
                if !strong.insert(h.clone()) {
                    return Err(format!(
                        "internal error: repeated strong generator at i={} j={} base={:?} strong_len={} h={:?}",
                        failure_level.0,
                        failure_level.1,
                        base,
                        strong.len(),
                        h
                    ));
                }
            }
            None => {
                let mut order: u128 = 1;
                for level in &levels {
                    order = order
                        .checked_mul(level.transversal.len() as u128)
                        .ok_or_else(|| {
                            "group order is too large for 128-bit arithmetic".to_string()
                        })?;
                }
                return Ok((base, strong.into_iter().collect(), levels, order));
            }
        }
    }
}

/// A finite permutation group generated by permutations of `{0, ..., n-1}`.
///
/// The group is represented by a base and strong generating set computed with
/// the Schreier-Sims algorithm. Order, membership, orbits, and stabilizers are
/// answered from that structure without enumerating the group elements.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct PermGroup {
    n: usize,
    generators: Vec<Perm>,
    base: Vec<usize>,
    strong: Vec<Images>,
    levels: Vec<StabilizerLevel>,
    order: u128,
}

impl PermGroup {
    pub(crate) fn from_generators(n: usize, generators: Vec<Perm>) -> PyResult<Self> {
        for g in &generators {
            if g.n() != n {
                return Err(PyValueError::new_err("generator has different n"));
            }
        }
        let (base, strong, levels, order) =
            build_perm_group(n, &generators).map_err(PyValueError::new_err)?;
        Ok(Self {
            n,
            generators,
            base,
            strong,
            levels,
            order,
        })
    }

    fn contains_images(&self, images: &[usize]) -> bool {
        if images.len() != self.n {
            return false;
        }
        let id = identity_images(self.n);
        let mut current = images.to_vec();
        for level in &self.levels {
            match level.transversal.get(&current[level.point]) {
                Some(u) => current = compose_images(&inverse_images(u), &current),
                None => return false,
            }
        }
        current == id
    }

    fn orbit_of(&self, point: usize) -> Vec<usize> {
        orbit_transversal(self.n, point, &self.strong).0
    }

    fn validate_point(&self, point: usize) -> PyResult<()> {
        if point >= self.n {
            return Err(PyValueError::new_err("point out of range"));
        }
        Ok(())
    }
}

#[pymethods]
impl PermGroup {
    #[new]
    /// Construct the group generated by `generators` acting on `{0, ..., n-1}`.
    pub fn new(n: usize, generators: Vec<Perm>) -> PyResult<Self> {
        Self::from_generators(n, generators)
    }

    /// Return the permutation degree.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Return the generators the group was constructed from.
    pub fn generators(&self) -> Vec<Perm> {
        self.generators.clone()
    }

    /// Return the Schreier-Sims base.
    pub fn base(&self) -> Vec<usize> {
        self.base.clone()
    }

    /// Return the strong generating set relative to the base.
    pub fn strong_generators(&self) -> Vec<Perm> {
        self.strong
            .iter()
            .map(|images| Perm::from_images_unchecked(self.n, images.clone()))
            .collect()
    }

    /// Return the exact group order. No elements are enumerated.
    pub fn order(&self) -> u128 {
        self.order
    }

    /// Return whether `perm` lies in the group.
    pub fn contains(&self, perm: &Perm) -> PyResult<bool> {
        if perm.n() != self.n {
            return Err(PyValueError::new_err("different n"));
        }
        Ok(self.contains_images(perm.images_ref()))
    }

    /// Return the orbit of `point` under the group.
    pub fn orbit(&self, point: usize) -> PyResult<Vec<usize>> {
        self.validate_point(point)?;
        Ok(self.orbit_of(point))
    }

    /// Return the orbits of the group on `{0, ..., n-1}`.
    pub fn orbits(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.n];
        let mut out: Vec<Vec<usize>> = Vec::new();
        for start in 0..self.n {
            if visited[start] {
                continue;
            }
            let orbit = self.orbit_of(start);
            for &x in &orbit {
                visited[x] = true;
            }
            out.push(orbit);
        }
        out
    }

    /// Return whether the group action is transitive.
    pub fn is_transitive(&self) -> bool {
        self.n == 0 || self.orbit_of(0).len() == self.n
    }

    /// Return the stabilizer of `point` as a generated subgroup.
    pub fn stabilizer(&self, point: usize) -> PyResult<PermGroup> {
        self.validate_point(point)?;
        let (orbit, transversal) = orbit_transversal(self.n, point, &self.strong);
        let id = identity_images(self.n);
        let mut generators: BTreeSet<Images> = BTreeSet::new();

        for x in &orbit {
            let ux = &transversal[x];
            for s in &self.strong {
                let y = s[*x];
                let uy = &transversal[&y];
                let schreier = compose_images(&compose_images(&inverse_images(uy), s), ux);
                if schreier != id {
                    generators.insert(schreier);
                }
            }
        }

        let generators: Vec<Perm> = generators
            .into_iter()
            .map(|images| Perm::from_images_unchecked(self.n, images))
            .collect();
        PermGroup::from_generators(self.n, generators)
    }

    /// Return the size of the stabilizer of `point` via orbit-stabilizer.
    pub fn stabilizer_size(&self, point: usize) -> PyResult<u128> {
        self.validate_point(point)?;
        Ok(self.order / self.orbit_of(point).len() as u128)
    }

    /// Return whether the group is abelian.
    pub fn is_abelian(&self) -> bool {
        for i in 0..self.generators.len() {
            for j in (i + 1)..self.generators.len() {
                let a = self.generators[i].images_ref();
                let b = self.generators[j].images_ref();
                if compose_images(a, b) != compose_images(b, a) {
                    return false;
                }
            }
        }
        true
    }

    /// Return the highest power of the prime `p` dividing the group order.
    pub fn p_part_order(&self, p: u64) -> PyResult<u128> {
        if !is_prime_u64(p) {
            return Err(PyValueError::new_err("p must be prime"));
        }
        let p = u128::from(p);
        let mut remaining = self.order;
        let mut part: u128 = 1;
        while remaining.is_multiple_of(p) {
            part *= p;
            remaining /= p;
        }
        Ok(part)
    }

    /// Return whether the group order is a power of the prime `p`.
    pub fn is_p_group(&self, p: u64) -> PyResult<bool> {
        Ok(self.p_part_order(p)? == self.order)
    }

    /// Enumerate the group elements as an explicit `PermSubgroup`, subject to `max_size`.
    #[pyo3(signature = (max_size=None))]
    pub fn as_explicit(&self, max_size: Option<usize>) -> PyResult<PermSubgroup> {
        let max = max_size.unwrap_or(4096);
        if self.order > max as u128 {
            return Err(PyValueError::new_err(format!(
                "group order {} exceeds max_size {}; pass a larger max_size to enumerate explicitly",
                self.order, max
            )));
        }
        PermSubgroup::from_generated(self.n, self.generators.clone(), max)
    }

    pub fn __len__(&self) -> PyResult<usize> {
        usize::try_from(self.order)
            .map_err(|_| PyValueError::new_err("group order does not fit in usize"))
    }

    pub fn __repr__(&self) -> String {
        format!("PermGroup(n={}, order={})", self.n, self.order)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(n: usize, cycles: &[&[usize]]) -> PermGroup {
        let generators = cycles
            .iter()
            .map(|cycle| Perm::cycle(n, cycle.to_vec()).unwrap())
            .collect();
        PermGroup::from_generators(n, generators).unwrap()
    }

    #[test]
    fn s4_has_order_24_without_enumeration() {
        let g = group(4, &[&[0, 1], &[0, 1, 2, 3]]);
        assert_eq!(g.order(), 24);
        assert!(g.is_transitive());
        assert!(!g.is_abelian());
        assert_eq!(g.stabilizer_size(0).unwrap(), 6);
    }

    #[test]
    fn membership_matches_named_elements() {
        let g = group(4, &[&[0, 1], &[0, 1, 2, 3]]);
        assert!(g.contains(&Perm::identity(4)).unwrap());
        assert!(
            g.contains(&Perm::new(4, vec![1, 0, 3, 2]).unwrap())
                .unwrap()
        );

        let c4 = group(4, &[&[0, 1, 2, 3]]);
        assert!(
            c4.contains(&Perm::cycle(4, vec![0, 1, 2, 3]).unwrap())
                .unwrap()
        );
        assert!(!c4.contains(&Perm::cycle(4, vec![0, 1]).unwrap()).unwrap());
    }

    #[test]
    fn cyclic_group_orbits_and_stabilizers() {
        let g = group(6, &[&[0, 1, 2, 3, 4, 5]]);
        assert_eq!(g.order(), 6);
        assert_eq!(g.orbit(0).unwrap(), vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(g.stabilizer_size(0).unwrap(), 1);
        assert!(g.is_abelian());
        assert_eq!(g.p_part_order(2).unwrap(), 2);
        assert_eq!(g.p_part_order(3).unwrap(), 3);
        assert!(!g.is_p_group(2).unwrap());
    }

    #[test]
    fn dihedral_group_detects_order_and_stabilizer() {
        let g = group(4, &[&[0, 1, 2, 3], &[1, 3]]);
        assert_eq!(g.order(), 8);
        let reflection = Perm::new(4, vec![0, 3, 2, 1]).unwrap();
        assert!(g.contains(&reflection).unwrap());
        assert_eq!(g.stabilizer_size(0).unwrap(), 2);
        assert_eq!(g.stabilizer(0).unwrap().order(), 2);
    }

    #[test]
    fn v4_is_abelian_and_transitive() {
        let g = PermGroup::from_generators(
            4,
            vec![
                Perm::new(4, vec![1, 0, 3, 2]).unwrap(),
                Perm::new(4, vec![2, 3, 0, 1]).unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(g.order(), 4);
        assert!(g.is_abelian());
        assert!(g.is_transitive());
    }

    #[test]
    fn large_symmetric_group_order_is_exact() {
        let n_cycle: Vec<usize> = (0..12).collect();
        let g = PermGroup::from_generators(
            12,
            vec![
                Perm::cycle(12, vec![0, 1]).unwrap(),
                Perm::cycle(12, n_cycle).unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(g.order(), 479_001_600);
        assert!(g.is_transitive());
        assert_eq!(g.stabilizer_size(0).unwrap(), 39_916_800);
    }

    #[test]
    fn explicit_conversion_respects_bounds() {
        let g = group(4, &[&[0, 1], &[0, 1, 2, 3]]);
        let explicit = g.as_explicit(Some(24)).unwrap();
        assert_eq!(explicit.order(), 24);
        assert!(g.as_explicit(Some(10)).is_err());
    }

    #[test]
    fn stabilizer_chain_order_matches_direct_enumeration() {
        let g = group(5, &[&[0, 1, 2], &[3, 4]]);
        let explicit = g.as_explicit(Some(100)).unwrap();
        assert_eq!(g.order() as usize, explicit.order());
        assert_eq!(g.orbits(), explicit.orbits());
        assert_eq!(
            g.stabilizer_size(0).unwrap() as usize,
            explicit.stabilizer(0).unwrap().order()
        );
    }
}

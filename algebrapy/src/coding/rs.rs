use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::field::fq::{Fq, FqElem};
use crate::field::poly_fp::PolyFp;

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct ReedSolomonCode {
    p: u64,
    n: u64,
    dimension: usize,
    primitive_modulus: PolyFp,
    generator_coeffs: Vec<Vec<u64>>,
}

#[pymethods]
impl ReedSolomonCode {
    /// Primitive narrow-sense Reed-Solomon code over GF(p^m) with length n=q-1.
    #[new]
    pub fn new(p: u64, primitive_modulus_coeffs: Vec<i128>, dimension: usize) -> PyResult<Self> {
        let primitive_modulus = PolyFp::new(p, primitive_modulus_coeffs)?;
        if !primitive_modulus.is_irreducible()? {
            return Err(PyValueError::new_err(
                "primitive modulus must be irreducible",
            ));
        }

        let fq = Fq::new(
            p,
            primitive_modulus
                .coeffs()
                .iter()
                .map(|&c| c as i128)
                .collect(),
        )?;
        let n = fq.size() - 1;
        if n == 0 {
            return Err(PyValueError::new_err(
                "field is too small for Reed-Solomon codes",
            ));
        }
        if dimension == 0 || dimension >= n as usize {
            return Err(PyValueError::new_err("dimension must satisfy 1 <= k < q-1"));
        }

        let alpha = fq.elem(vec![0, 1])?;
        if fq.mul_order(&alpha)? != n {
            return Err(PyValueError::new_err(
                "primitive modulus does not make x a primitive element of GF(p^m)",
            ));
        }

        let generator = build_rs_generator(&fq, &alpha, n - dimension as u64)?;
        Ok(Self {
            p,
            n,
            dimension,
            primitive_modulus,
            generator_coeffs: generator.iter().map(FqElem::coeffs).collect(),
        })
    }

    pub fn length(&self) -> u64 {
        self.n
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn designed_distance(&self) -> u64 {
        self.n - self.dimension as u64 + 1
    }

    pub fn correction_capacity(&self) -> usize {
        ((self.n as usize) - self.dimension) / 2
    }

    pub fn field(&self) -> PyResult<Fq> {
        self.make_field()
    }

    pub fn generator_poly(&self) -> PyResult<Vec<FqElem>> {
        self.make_generator_poly()
    }

    pub fn encode_systematic(&self, message: Vec<FqElem>) -> PyResult<Vec<FqElem>> {
        self.check_message_len(&message)?;
        let fq = self.make_field()?;
        let generator = self.make_generator_poly()?;
        check_symbols(&fq, &message)?;

        let shifted = fq_poly_shift(&fq, &message, generator.len() - 1);
        let remainder = fq_poly_mod(&fq, &shifted, &generator)?;
        let codeword = fq_poly_sub(&fq, &shifted, &remainder)?;
        Ok(fq_poly_pad(codeword, self.n as usize, &fq.zero()))
    }

    pub fn is_codeword(&self, word: Vec<FqElem>) -> PyResult<bool> {
        self.check_word_len(&word)?;
        let fq = self.make_field()?;
        let generator = self.make_generator_poly()?;
        check_symbols(&fq, &word)?;
        Ok(fq_poly_mod(&fq, &word, &generator)?
            .iter()
            .all(FqElem::is_zero))
    }

    pub fn syndromes(&self, word: Vec<FqElem>) -> PyResult<Vec<FqElem>> {
        self.check_word_len(&word)?;
        let fq = self.make_field()?;
        let alpha = self.primitive_element(&fq)?;
        check_symbols(&fq, &word)?;
        let parity_symbols = (self.n as usize) - self.dimension;
        let mut out = Vec::with_capacity(parity_symbols);
        for i in 1..=parity_symbols {
            let point = fq.pow(&alpha, i as i128)?;
            out.push(fq_poly_eval(&fq, &word, &point)?);
        }
        Ok(out)
    }

    pub fn decode(&self, word: Vec<FqElem>) -> PyResult<Vec<FqElem>> {
        self.check_word_len(&word)?;
        let fq = self.make_field()?;
        let alpha = self.primitive_element(&fq)?;
        check_symbols(&fq, &word)?;

        let syndromes = self.syndromes(word.clone())?;
        if syndromes.iter().all(FqElem::is_zero) {
            return Ok(word);
        }

        let sigma = berlekamp_massey(&fq, &syndromes)?;
        let locator_degree = fq_poly_degree(&sigma);
        if locator_degree == 0 || locator_degree > self.correction_capacity() {
            return Err(PyValueError::new_err(
                "received word appears to have too many symbol errors",
            ));
        }

        let error_locations = rs_chien_search(&fq, &alpha, self.n, &sigma)?;
        if error_locations.len() != locator_degree {
            return Err(PyValueError::new_err(
                "failed to locate the expected number of Reed-Solomon errors",
            ));
        }

        let magnitudes = solve_error_magnitudes(&fq, &alpha, &syndromes, &error_locations)?;
        let mut corrected = word;
        for ((pos, _root), magnitude) in error_locations.into_iter().zip(magnitudes.into_iter()) {
            corrected[pos as usize] = fq.sub(&corrected[pos as usize], &magnitude)?;
        }

        if !self.is_codeword(corrected.clone())? {
            return Err(PyValueError::new_err(
                "Reed-Solomon decoding failed verification",
            ));
        }
        Ok(corrected)
    }

    pub fn extract_systematic_message(&self, codeword: Vec<FqElem>) -> PyResult<Vec<FqElem>> {
        self.check_word_len(&codeword)?;
        if !self.is_codeword(codeword.clone())? {
            return Err(PyValueError::new_err("word is not a Reed-Solomon codeword"));
        }
        let parity_symbols = (self.n as usize) - self.dimension;
        let zero = self.make_field()?.zero();
        let mut out = if codeword.len() <= parity_symbols {
            vec![]
        } else {
            codeword[parity_symbols..].to_vec()
        };
        out.resize(self.dimension, zero);
        Ok(out)
    }

    pub fn decode_message(&self, word: Vec<FqElem>) -> PyResult<Vec<FqElem>> {
        let corrected = self.decode(word)?;
        self.extract_systematic_message(corrected)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "ReedSolomonCode(p={}, n={}, k={}, d={})",
            self.p,
            self.n,
            self.dimension,
            self.designed_distance()
        )
    }
}

impl ReedSolomonCode {
    fn make_field(&self) -> PyResult<Fq> {
        Fq::new(
            self.p,
            self.primitive_modulus
                .coeffs()
                .iter()
                .map(|&c| c as i128)
                .collect(),
        )
    }

    fn primitive_element(&self, fq: &Fq) -> PyResult<FqElem> {
        fq.elem(vec![0, 1])
    }

    fn make_generator_poly(&self) -> PyResult<Vec<FqElem>> {
        let fq = self.make_field()?;
        self.generator_coeffs
            .iter()
            .map(|c| fq.elem(c.iter().map(|&x| x as i128).collect()))
            .collect()
    }

    fn check_message_len(&self, message: &[FqElem]) -> PyResult<()> {
        if message.len() > self.dimension {
            Err(PyValueError::new_err(
                "message is too long for the Reed-Solomon code dimension",
            ))
        } else {
            Ok(())
        }
    }

    fn check_word_len(&self, word: &[FqElem]) -> PyResult<()> {
        if word.len() > self.n as usize {
            Err(PyValueError::new_err(
                "word is longer than the Reed-Solomon code length",
            ))
        } else {
            Ok(())
        }
    }
}

fn build_rs_generator(fq: &Fq, alpha: &FqElem, parity_symbols: u64) -> PyResult<Vec<FqElem>> {
    let mut g = vec![fq.one()];
    for i in 1..=parity_symbols {
        let root = fq.pow(alpha, i as i128)?;
        let factor = vec![fq.neg(&root)?, fq.one()];
        g = fq_poly_mul(fq, &g, &factor)?;
    }
    Ok(fq_poly_trim(g))
}

fn check_symbols(fq: &Fq, symbols: &[FqElem]) -> PyResult<()> {
    for s in symbols {
        fq.add(s, &fq.zero())?;
    }
    Ok(())
}

fn fq_poly_trim(mut coeffs: Vec<FqElem>) -> Vec<FqElem> {
    while coeffs.last().is_some_and(FqElem::is_zero) {
        coeffs.pop();
    }
    coeffs
}

fn fq_poly_degree(coeffs: &[FqElem]) -> usize {
    coeffs.iter().rposition(|c| !c.is_zero()).unwrap_or(0)
}

fn fq_poly_shift(fq: &Fq, coeffs: &[FqElem], shift: usize) -> Vec<FqElem> {
    if coeffs.is_empty() {
        return vec![];
    }
    let mut out = vec![fq.zero(); shift];
    out.extend_from_slice(coeffs);
    out
}

fn fq_poly_pad(mut coeffs: Vec<FqElem>, len: usize, zero: &FqElem) -> Vec<FqElem> {
    while coeffs.len() < len {
        coeffs.push(zero.clone());
    }
    coeffs
}

fn fq_poly_get<'a>(coeffs: &'a [FqElem], idx: usize, zero: &'a FqElem) -> &'a FqElem {
    coeffs.get(idx).unwrap_or(zero)
}

fn fq_poly_sub(fq: &Fq, a: &[FqElem], b: &[FqElem]) -> PyResult<Vec<FqElem>> {
    let n = a.len().max(b.len());
    let zero = fq.zero();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(fq.sub(fq_poly_get(a, i, &zero), fq_poly_get(b, i, &zero))?);
    }
    Ok(fq_poly_trim(out))
}

fn fq_poly_mul(fq: &Fq, a: &[FqElem], b: &[FqElem]) -> PyResult<Vec<FqElem>> {
    if a.is_empty() || b.is_empty() {
        return Ok(vec![]);
    }
    let zero = fq.zero();
    let mut out = vec![zero.clone(); a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            let term = fq.mul(ai, bj)?;
            out[i + j] = fq.add(&out[i + j], &term)?;
        }
    }
    Ok(fq_poly_trim(out))
}

fn fq_poly_scale_shift(
    fq: &Fq,
    poly: &[FqElem],
    scale: &FqElem,
    shift: usize,
) -> PyResult<Vec<FqElem>> {
    let mut out = vec![fq.zero(); poly.len() + shift];
    for (i, coeff) in poly.iter().enumerate() {
        out[i + shift] = fq.mul(coeff, scale)?;
    }
    Ok(fq_poly_trim(out))
}

fn fq_poly_div_rem(
    fq: &Fq,
    numer: &[FqElem],
    denom: &[FqElem],
) -> PyResult<(Vec<FqElem>, Vec<FqElem>)> {
    let mut r = fq_poly_trim(numer.to_vec());
    let d = fq_poly_trim(denom.to_vec());
    if d.is_empty() {
        return Err(PyValueError::new_err("division by zero polynomial"));
    }
    if r.len() < d.len() {
        return Ok((vec![], r));
    }
    let zero = fq.zero();
    let mut q = vec![zero.clone(); r.len() - d.len() + 1];
    let d_lc_inv = fq.inv(d.last().unwrap())?;

    while !r.is_empty() && r.len() >= d.len() {
        let shift = r.len() - d.len();
        let factor = fq.mul(r.last().unwrap(), &d_lc_inv)?;
        q[shift] = factor.clone();
        for (i, dc) in d.iter().enumerate() {
            let idx = i + shift;
            let term = fq.mul(&factor, dc)?;
            r[idx] = fq.sub(&r[idx], &term)?;
        }
        r = fq_poly_trim(r);
    }
    Ok((fq_poly_trim(q), r))
}

fn fq_poly_mod(fq: &Fq, numer: &[FqElem], denom: &[FqElem]) -> PyResult<Vec<FqElem>> {
    Ok(fq_poly_div_rem(fq, numer, denom)?.1)
}

fn fq_poly_eval(fq: &Fq, coeffs: &[FqElem], x: &FqElem) -> PyResult<FqElem> {
    let mut acc = fq.zero();
    for coeff in coeffs.iter().rev() {
        acc = fq.mul(&acc, x)?;
        acc = fq.add(&acc, coeff)?;
    }
    Ok(acc)
}

fn berlekamp_massey(fq: &Fq, syndromes: &[FqElem]) -> PyResult<Vec<FqElem>> {
    let zero = fq.zero();
    let one = fq.one();
    let mut c = vec![one.clone()];
    let mut b = vec![one];
    let mut l = 0usize;
    let mut m = 1usize;
    let mut bb = fq.one();

    for n in 0..syndromes.len() {
        let mut d = syndromes[n].clone();
        for i in 1..=l {
            let term = fq.mul(fq_poly_get(&c, i, &zero), &syndromes[n - i])?;
            d = fq.add(&d, &term)?;
        }
        if d.is_zero() {
            m += 1;
            continue;
        }
        let factor = fq.mul(&d, &fq.inv(&bb)?)?;
        let correction = fq_poly_scale_shift(fq, &b, &factor, m)?;
        let t = fq_poly_sub(fq, &c, &correction)?;
        if 2 * l <= n {
            l = n + 1 - l;
            b = c;
            bb = d;
            c = t;
            m = 1;
        } else {
            c = t;
            m += 1;
        }
    }
    Ok(fq_poly_trim(c))
}

fn rs_chien_search(
    fq: &Fq,
    alpha: &FqElem,
    n: u64,
    sigma: &[FqElem],
) -> PyResult<Vec<(u64, FqElem)>> {
    let mut out = Vec::new();
    for j in 0..n {
        let root = fq.pow(alpha, -(j as i128))?;
        if fq_poly_eval(fq, sigma, &root)?.is_zero() {
            out.push((j, root));
        }
    }
    Ok(out)
}

fn solve_error_magnitudes(
    fq: &Fq,
    alpha: &FqElem,
    syndromes: &[FqElem],
    error_locations: &[(u64, FqElem)],
) -> PyResult<Vec<FqElem>> {
    let nu = error_locations.len();
    if nu == 0 {
        return Ok(vec![]);
    }
    let mut matrix = vec![vec![fq.zero(); nu + 1]; nu];
    for row in 0..nu {
        for (col, (pos, _root)) in error_locations.iter().enumerate() {
            let x = fq.pow(alpha, *pos as i128)?;
            matrix[row][col] = fq.pow(&x, (row + 1) as i128)?;
        }
        matrix[row][nu] = syndromes[row].clone();
    }
    gaussian_elimination(fq, matrix)
}

fn gaussian_elimination(fq: &Fq, mut a: Vec<Vec<FqElem>>) -> PyResult<Vec<FqElem>> {
    let n = a.len();
    for col in 0..n {
        let pivot = (col..n)
            .find(|&row| !a[row][col].is_zero())
            .ok_or_else(|| PyValueError::new_err("singular Reed-Solomon decoding system"))?;
        if pivot != col {
            a.swap(pivot, col);
        }

        let inv = fq.inv(&a[col][col])?;
        for j in col..=n {
            a[col][j] = fq.mul(&a[col][j], &inv)?;
        }

        for row in 0..n {
            if row == col || a[row][col].is_zero() {
                continue;
            }
            let factor = a[row][col].clone();
            for j in col..=n {
                let term = fq.mul(&factor, &a[col][j])?;
                a[row][j] = fq.sub(&a[row][j], &term)?;
            }
        }
    }
    Ok((0..n).map(|i| a[i][n].clone()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gf8() -> ReedSolomonCode {
        ReedSolomonCode::new(2, vec![1, 1, 0, 1], 3).unwrap()
    }

    #[test]
    fn rs_parameters_and_generator_degree() {
        let rs = gf8();
        assert_eq!(rs.length(), 7);
        assert_eq!(rs.dimension(), 3);
        assert_eq!(rs.designed_distance(), 5);
        assert_eq!(rs.correction_capacity(), 2);
        assert_eq!(rs.generator_poly().unwrap().len(), 5);
    }

    #[test]
    fn rs_systematic_round_trip() {
        let rs = gf8();
        let fq = rs.field().unwrap();
        let msg = vec![
            fq.elem(vec![1]).unwrap(),
            fq.elem(vec![0, 1]).unwrap(),
            fq.elem(vec![1, 1]).unwrap(),
        ];
        let codeword = rs.encode_systematic(msg.clone()).unwrap();
        assert!(rs.is_codeword(codeword.clone()).unwrap());
        let recovered = rs.extract_systematic_message(codeword).unwrap();
        assert_eq!(
            recovered.iter().map(FqElem::coeffs).collect::<Vec<_>>(),
            msg.iter().map(FqElem::coeffs).collect::<Vec<_>>()
        );
    }

    #[test]
    fn rs_decodes_two_symbol_errors() {
        let rs = gf8();
        let fq = rs.field().unwrap();
        let msg = vec![
            fq.elem(vec![1]).unwrap(),
            fq.elem(vec![0, 1]).unwrap(),
            fq.elem(vec![1, 1]).unwrap(),
        ];
        let mut received = rs.encode_systematic(msg.clone()).unwrap();
        received[1] = fq.add(&received[1], &fq.elem(vec![1]).unwrap()).unwrap();
        received[5] = fq.add(&received[5], &fq.elem(vec![1, 1]).unwrap()).unwrap();
        let decoded = rs.decode(received).unwrap();
        let recovered = rs.extract_systematic_message(decoded).unwrap();
        assert_eq!(
            recovered.iter().map(FqElem::coeffs).collect::<Vec<_>>(),
            msg.iter().map(FqElem::coeffs).collect::<Vec<_>>()
        );
    }

    #[test]
    fn rs_rejects_too_many_errors() {
        let rs = gf8();
        let fq = rs.field().unwrap();
        let msg = vec![
            fq.elem(vec![1]).unwrap(),
            fq.elem(vec![0, 1]).unwrap(),
            fq.elem(vec![1, 1]).unwrap(),
        ];
        let mut received = rs.encode_systematic(msg).unwrap();
        received[0] = fq.add(&received[0], &fq.elem(vec![1]).unwrap()).unwrap();
        received[1] = fq.add(&received[1], &fq.elem(vec![1, 1]).unwrap()).unwrap();
        received[2] = fq.add(&received[2], &fq.elem(vec![0, 1]).unwrap()).unwrap();
        assert!(rs.decode(received).is_err());
    }

    #[test]
    fn rs_syndromes_vanish_on_codeword() {
        let rs = gf8();
        let fq = rs.field().unwrap();
        let msg = vec![fq.elem(vec![1]).unwrap(), fq.zero(), fq.one()];
        let codeword = rs.encode_systematic(msg).unwrap();
        assert!(rs.syndromes(codeword).unwrap().iter().all(FqElem::is_zero));
    }

    #[test]
    fn rs_decodes_all_two_symbol_errors_for_representative_messages() {
        let rs = gf8();
        let fq = rs.field().unwrap();
        let elements = fq.elements(Some(16)).unwrap();
        let n = rs.length() as usize;
        let messages = vec![
            vec![fq.zero(), fq.zero(), fq.zero()],
            vec![fq.one(), fq.zero(), fq.zero()],
            vec![fq.zero(), fq.one(), fq.zero()],
            vec![fq.zero(), fq.zero(), fq.one()],
            vec![fq.elem(vec![0, 1]).unwrap(), fq.zero(), fq.one()],
            vec![
                fq.elem(vec![1, 1]).unwrap(),
                fq.elem(vec![1, 0, 1]).unwrap(),
                fq.zero(),
            ],
            vec![
                fq.elem(vec![1]).unwrap(),
                fq.elem(vec![0, 1]).unwrap(),
                fq.elem(vec![1, 1]).unwrap(),
            ],
            vec![
                fq.elem(vec![1, 1, 1]).unwrap(),
                fq.elem(vec![1, 0]).unwrap(),
                fq.elem(vec![0, 0, 1]).unwrap(),
            ],
        ];

        for msg in messages {
            let codeword = rs.encode_systematic(msg.clone()).unwrap();

            for i in 0..n {
                for j in (i + 1)..n {
                    for e1 in &elements {
                        if e1.is_zero() {
                            continue;
                        }
                        for e2 in &elements {
                            if e2.is_zero() {
                                continue;
                            }
                            let mut received = codeword.clone();
                            received[i] = fq.add(&received[i], e1).unwrap();
                            received[j] = fq.add(&received[j], e2).unwrap();
                            let decoded = rs.decode(received).unwrap();
                            let recovered = rs.extract_systematic_message(decoded).unwrap();
                            assert_eq!(
                                recovered.iter().map(FqElem::coeffs).collect::<Vec<_>>(),
                                msg.iter().map(FqElem::coeffs).collect::<Vec<_>>()
                            );
                        }
                    }
                }
            }
        }
    }
}

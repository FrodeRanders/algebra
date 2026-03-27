use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::{BTreeSet, HashSet};

use crate::field::fq::{Fq, FqElem};
use crate::field::poly_fp::PolyFp;

#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct BinaryBchCode {
    m: u32,
    n: u64,
    designed_distance: u64,
    primitive_modulus: PolyFp,
    generator: PolyFp,
}

#[pymethods]
impl BinaryBchCode {
    /// Primitive narrow-sense binary BCH code of length n = 2^m - 1.
    ///
    /// `primitive_modulus_coeffs` defines the irreducible polynomial for GF(2^m)
    /// in low-to-high order. The code uses alpha = x mod f(x) as a primitive
    /// element candidate and validates that it has multiplicative order n.
    #[new]
    pub fn new(m: u32, primitive_modulus_coeffs: Vec<i128>, designed_distance: u64) -> PyResult<Self> {
        if m == 0 {
            return Err(PyValueError::new_err("m must be >= 1"));
        }
        if m >= 63 {
            return Err(PyValueError::new_err("m must be < 63"));
        }
        let n = (1u64 << m) - 1;
        if designed_distance < 2 {
            return Err(PyValueError::new_err("designed_distance must be >= 2"));
        }
        if designed_distance > n {
            return Err(PyValueError::new_err(
                "designed_distance must be <= code length 2^m - 1",
            ));
        }

        let primitive_modulus = PolyFp::new(2, primitive_modulus_coeffs)?;
        if primitive_modulus.degree() != m as i64 {
            return Err(PyValueError::new_err(
                "primitive modulus degree must equal m",
            ));
        }
        if !primitive_modulus.is_irreducible()? {
            return Err(PyValueError::new_err(
                "primitive modulus must be irreducible over F2",
            ));
        }

        let fq = Fq::new(
            2,
            primitive_modulus
                .coeffs()
                .iter()
                .map(|&c| c as i128)
                .collect(),
        )?;
        let alpha = fq.elem(vec![0, 1])?;
        if fq.mul_order(&alpha)? != n {
            return Err(PyValueError::new_err(
                "primitive modulus does not make x a primitive element of GF(2^m)",
            ));
        }

        let generator = build_binary_bch_generator(&fq, &alpha, n, designed_distance)?;

        Ok(Self {
            m,
            n,
            designed_distance,
            primitive_modulus,
            generator,
        })
    }

    pub fn length(&self) -> u64 {
        self.n
    }

    pub fn designed_distance(&self) -> u64 {
        self.designed_distance
    }

    pub fn generator_poly(&self) -> PolyFp {
        self.generator.clone()
    }

    pub fn generator_degree(&self) -> usize {
        self.generator.coeffs().len().saturating_sub(1)
    }

    pub fn dimension(&self) -> usize {
        (self.n as usize) - self.generator_degree()
    }

    pub fn correction_capacity(&self) -> u64 {
        (self.designed_distance - 1) / 2
    }

    /// Encode a binary message polynomial m(x) as c(x) = m(x) g(x).
    /// Coefficients are given low-to-high and must be 0/1.
    pub fn encode(&self, message_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        if message_bits.len() > self.dimension() {
            return Err(PyValueError::new_err(
                "message is too long for the BCH code dimension",
            ));
        }
        let msg = binary_poly_from_bits(&message_bits)?;
        let codeword = msg.mul(&self.generator)?;
        let bits = bits_from_binary_poly(&codeword)?;
        if bits.len() > self.n as usize {
            return Err(PyValueError::new_err(
                "encoded polynomial exceeds the code length",
            ));
        }
        Ok(pad_bits(bits, self.n as usize))
    }

    /// Systematic cyclic encoding c(x) = x^(n-k)m(x) + r(x), where
    /// r(x) = x^(n-k)m(x) mod g(x).
    pub fn encode_systematic(&self, message_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        if message_bits.len() > self.dimension() {
            return Err(PyValueError::new_err(
                "message is too long for the BCH code dimension",
            ));
        }
        let msg = binary_poly_from_bits(&message_bits)?;
        let shifted = shift_poly(&msg, self.generator_degree())?;
        let parity = shifted.modulo(&self.generator)?;
        let codeword = shifted.add(&parity)?;
        let bits = bits_from_binary_poly(&codeword)?;
        if bits.len() > self.n as usize {
            return Err(PyValueError::new_err(
                "systematic encoding exceeded the BCH code length",
            ));
        }
        Ok(pad_bits(bits, self.n as usize))
    }

    /// Recover the message polynomial from a valid codeword by dividing by g(x).
    pub fn extract_message(&self, codeword_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        self.check_word_length(&codeword_bits)?;
        let word = binary_poly_from_bits(&codeword_bits)?;
        let (quotient, remainder) = word.div_rem(&self.generator)?;
        if !remainder.is_zero() {
            return Err(PyValueError::new_err(
                "word is not divisible by the BCH generator polynomial",
            ));
        }
        let bits = bits_from_binary_poly(&quotient)?;
        if bits.len() > self.dimension() {
            return Err(PyValueError::new_err(
                "recovered message exceeds the BCH code dimension",
            ));
        }
        Ok(pad_bits(bits, self.dimension()))
    }

    /// Recover the message bits from a valid systematic codeword by reading the
    /// high-degree message block after validating the codeword.
    pub fn extract_systematic_message(&self, codeword_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        self.check_word_length(&codeword_bits)?;
        if !self.is_codeword(codeword_bits.clone())? {
            return Err(PyValueError::new_err(
                "word is not a valid BCH codeword",
            ));
        }
        let parity_len = self.generator_degree();
        let mut out = if codeword_bits.len() <= parity_len {
            vec![]
        } else {
            codeword_bits[parity_len..].to_vec()
        };
        out.resize(self.dimension(), 0);
        Ok(out)
    }

    /// Check the cyclic divisibility condition r(x) mod g(x) == 0.
    pub fn is_codeword(&self, word_bits: Vec<u64>) -> PyResult<bool> {
        self.check_word_length(&word_bits)?;
        let word = binary_poly_from_bits(&word_bits)?;
        Ok(word.modulo(&self.generator)?.is_zero())
    }

    /// Return the syndromes [r(alpha), r(alpha^2), ..., r(alpha^(d-1))].
    pub fn syndromes(&self, word_bits: Vec<u64>) -> PyResult<Vec<Vec<u64>>> {
        self.check_word_length(&word_bits)?;
        let fq = self.extension_field()?;
        let alpha = fq.elem(vec![0, 1])?;
        let word = word_bits;
        let mut out = Vec::with_capacity((self.designed_distance - 1) as usize);
        for i in 1..self.designed_distance {
            let root = fq.pow(&alpha, i as i128)?;
            let value = eval_binary_poly_at(&fq, &word, &root)?;
            out.push(value.coeffs());
        }
        Ok(out)
    }

    pub fn parity_check(&self, word_bits: Vec<u64>) -> PyResult<bool> {
        Ok(self.syndromes(word_bits)?.iter().all(|s| s.is_empty()))
    }

    /// Decode with a binary BCH syndrome decoder using Berlekamp-Massey and
    /// Chien search.
    ///
    /// The result is verified with a parity check, which catches many
    /// undecodable patterns. As with standard bounded-distance BCH decoding,
    /// words with more than t errors may still miscorrect to a different valid
    /// codeword; this method does not claim guaranteed detection beyond t.
    pub fn decode(&self, word_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        self.check_word_length(&word_bits)?;
        let t = self.correction_capacity() as usize;
        if t == 0 {
            return Ok(word_bits);
        }

        let fq = self.extension_field()?;
        let alpha = fq.elem(vec![0, 1])?;
        let syndromes = syndrome_elements(&fq, &alpha, &word_bits, 2 * t)?;
        if syndromes.iter().all(|s| s.is_zero()) {
            return Ok(word_bits);
        }

        let sigma = berlekamp_massey(&fq, &syndromes)?;
        let locator_degree = fq_poly_degree(&sigma);
        if locator_degree == 0 {
            return Err(PyValueError::new_err(
                "failed to build a nontrivial error locator polynomial",
            ));
        }
        if locator_degree > t {
            return Err(PyValueError::new_err(
                "received word appears to have more errors than the BCH decoder can correct",
            ));
        }

        let error_positions = chien_search(&fq, &alpha, self.n, &sigma)?;
        if error_positions.len() != locator_degree {
            return Err(PyValueError::new_err(
                "failed to locate the expected number of BCH error positions",
            ));
        }

        let mut corrected = word_bits;
        for pos in error_positions {
            corrected[pos as usize] ^= 1;
        }
        if !self.parity_check(corrected.clone())? {
            return Err(PyValueError::new_err(
                "BCH decoding failed verification; the received word may have too many errors",
            ));
        }
        Ok(corrected)
    }

    /// Decode a received word and then extract the underlying message
    /// polynomial from the corrected codeword.
    pub fn decode_message(&self, word_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        let corrected = self.decode(word_bits)?;
        self.extract_message(corrected)
    }

    pub fn decode_systematic_message(&self, word_bits: Vec<u64>) -> PyResult<Vec<u64>> {
        let corrected = self.decode(word_bits)?;
        self.extract_systematic_message(corrected)
    }

    /// Generator matrix formed by systematic encoding of basis messages.
    pub fn generator_matrix(&self) -> PyResult<Vec<Vec<u64>>> {
        let k = self.dimension();
        let mut rows = Vec::with_capacity(k);
        for i in 0..k {
            let mut msg = vec![0; k];
            msg[i] = 1;
            rows.push(self.encode_systematic(msg)?);
        }
        Ok(rows)
    }

    /// Parity-check matrix compatible with the systematic generator matrix
    /// returned by generator_matrix().
    pub fn parity_check_matrix(&self) -> PyResult<Vec<Vec<u64>>> {
        let g = self.generator_matrix()?;
        let r = self.generator_degree();
        let k = self.dimension();
        let mut h = vec![vec![0; self.n as usize]; r];
        for i in 0..r {
            h[i][i] = 1;
        }
        for row in 0..k {
            for col in 0..r {
                h[col][r + row] = g[row][col];
            }
        }
        Ok(h)
    }

    /// Parameters of the shortened systematic BCH code obtained by fixing the
    /// last `shorten_by` message coordinates to zero and puncturing those
    /// positions.
    pub fn shortened_parameters(&self, shorten_by: usize) -> PyResult<(u64, usize)> {
        let k = self.dimension();
        if shorten_by > k {
            return Err(PyValueError::new_err(
                "cannot shorten by more than the code dimension",
            ));
        }
        Ok((self.n - shorten_by as u64, k - shorten_by))
    }

    /// Shortened systematic encoding. The last `shorten_by` systematic message
    /// coordinates are fixed to zero and then punctured from the codeword.
    pub fn encode_shortened_systematic(
        &self,
        message_bits: Vec<u64>,
        shorten_by: usize,
    ) -> PyResult<Vec<u64>> {
        let (_n_short, k_short) = self.shortened_parameters(shorten_by)?;
        if message_bits.len() > k_short {
            return Err(PyValueError::new_err(
                "message is too long for the shortened BCH code dimension",
            ));
        }
        let mut full_message = message_bits;
        full_message.resize(self.dimension(), 0);
        let full_codeword = self.encode_systematic(full_message)?;
        Ok(full_codeword[..(self.n as usize - shorten_by)].to_vec())
    }

    /// Recover the message from a shortened systematic codeword by reattaching
    /// the punctured zero coordinates and using the full systematic extractor.
    pub fn extract_shortened_systematic_message(
        &self,
        codeword_bits: Vec<u64>,
        shorten_by: usize,
    ) -> PyResult<Vec<u64>> {
        let (n_short, k_short) = self.shortened_parameters(shorten_by)?;
        if codeword_bits.len() > n_short as usize {
            return Err(PyValueError::new_err(
                "word is longer than the shortened BCH code length",
            ));
        }
        let mut full_word = codeword_bits;
        full_word.resize(self.n as usize, 0);
        let mut msg = self.extract_systematic_message(full_word)?;
        msg.truncate(k_short);
        Ok(msg)
    }

    pub fn decode_shortened_systematic_message(
        &self,
        word_bits: Vec<u64>,
        shorten_by: usize,
    ) -> PyResult<Vec<u64>> {
        let (n_short, k_short) = self.shortened_parameters(shorten_by)?;
        if word_bits.len() > n_short as usize {
            return Err(PyValueError::new_err(
                "word is longer than the shortened BCH code length",
            ));
        }
        let mut full_word = word_bits;
        full_word.resize(self.n as usize, 0);
        let decoded = self.decode_systematic_message(full_word)?;
        Ok(decoded[..k_short].to_vec())
    }

    /// Enumerate all codewords when the dimension is small enough.
    #[pyo3(signature = (max_dimension=None))]
    pub fn codewords(&self, max_dimension: Option<usize>) -> PyResult<Vec<Vec<u64>>> {
        let limit = max_dimension.unwrap_or(16);
        let k = self.dimension();
        if k > limit {
            return Err(PyValueError::new_err(
                "code dimension too large to enumerate all codewords",
            ));
        }
        let count = 1usize << k;
        let mut out = Vec::with_capacity(count);
        for mask in 0..count {
            let msg = bits_from_mask(mask, k);
            out.push(self.encode_systematic(msg)?);
        }
        Ok(out)
    }

    /// Minimum nonzero Hamming weight, by exhaustive enumeration when k is small.
    #[pyo3(signature = (max_dimension=None))]
    pub fn minimum_distance(&self, max_dimension: Option<usize>) -> PyResult<usize> {
        let codewords = self.codewords(max_dimension)?;
        codewords
            .into_iter()
            .map(|cw| hamming_weight(&cw))
            .filter(|&w| w > 0)
            .min()
            .ok_or_else(|| PyValueError::new_err("failed to determine minimum distance"))
    }

    /// Weight distribution A_w for small codes, returned as counts indexed by weight.
    #[pyo3(signature = (max_dimension=None))]
    pub fn weight_distribution(&self, max_dimension: Option<usize>) -> PyResult<Vec<u64>> {
        let codewords = self.codewords(max_dimension)?;
        let mut counts = vec![0u64; (self.n as usize) + 1];
        for cw in codewords {
            counts[hamming_weight(&cw)] += 1;
        }
        Ok(counts)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "BinaryBchCode(m={}, n={}, d={}, k={}, g={:?})",
            self.m,
            self.n,
            self.designed_distance,
            self.dimension(),
            self.generator.coeffs()
        )
    }
}

impl BinaryBchCode {
    fn extension_field(&self) -> PyResult<Fq> {
        Fq::new(
            2,
            self.primitive_modulus
                .coeffs()
                .iter()
                .map(|&c| c as i128)
                .collect(),
        )
    }

    fn check_word_length(&self, word_bits: &[u64]) -> PyResult<()> {
        if word_bits.len() > self.n as usize {
            return Err(PyValueError::new_err(
                "word is longer than the BCH code length",
            ));
        }
        if word_bits.iter().any(|&b| b > 1) {
            return Err(PyValueError::new_err(
                "binary words must contain only 0/1 coefficients",
            ));
        }
        Ok(())
    }
}

fn build_binary_bch_generator(
    fq: &Fq,
    alpha: &FqElem,
    n: u64,
    designed_distance: u64,
) -> PyResult<PolyFp> {
    let mut covered: HashSet<u64> = HashSet::new();
    let mut generator = PolyFp::new(2, vec![1])?;

    for i in 1..designed_distance {
        let coset = cyclotomic_coset(i % n, n);
        let leader = *coset.first().unwrap();
        if covered.contains(&leader) {
            continue;
        }
        for &x in &coset {
            covered.insert(x);
        }
        let min_poly = minimal_polynomial_from_coset(fq, alpha, &coset)?;
        generator = generator.mul(&min_poly)?.monic()?;
    }

    Ok(generator.monic()?)
}

fn cyclotomic_coset(start: u64, n: u64) -> Vec<u64> {
    let mut seen = BTreeSet::new();
    let mut x = start % n;
    loop {
        if !seen.insert(x) {
            break;
        }
        x = (2 * x) % n;
    }
    seen.into_iter().collect()
}

fn minimal_polynomial_from_coset(fq: &Fq, alpha: &FqElem, coset: &[u64]) -> PyResult<PolyFp> {
    let mut coeffs: Vec<FqElem> = vec![fq.one()];

    for &exp in coset {
        let beta = fq.pow(alpha, exp as i128)?;
        let mut next = vec![fq.zero(); coeffs.len() + 1];
        for (i, coeff) in coeffs.iter().enumerate() {
            next[i + 1] = fq.add(&next[i + 1], coeff)?;
            let prod = fq.mul(coeff, &beta)?;
            next[i] = fq.sub(&next[i], &prod)?;
        }
        coeffs = next;
    }

    let mut binary_coeffs = Vec::with_capacity(coeffs.len());
    for coeff in coeffs {
        if coeff.is_zero() {
            binary_coeffs.push(0);
        } else if coeff.coeffs() == vec![1] {
            binary_coeffs.push(1);
        } else {
            return Err(PyValueError::new_err(
                "minimal polynomial did not land in F2; construction failed",
            ));
        }
    }

    PolyFp::new(2, binary_coeffs.iter().map(|&c| c as i128).collect())
}

fn eval_binary_poly_at(fq: &Fq, coeffs: &[u64], point: &FqElem) -> PyResult<FqElem> {
    let mut acc = fq.zero();
    let mut power = fq.one();

    for &bit in coeffs {
        if bit > 1 {
            return Err(PyValueError::new_err(
                "binary words must contain only 0/1 coefficients",
            ));
        }
        if bit == 1 {
            acc = fq.add(&acc, &power)?;
        }
        power = fq.mul(&power, point)?;
    }

    Ok(acc)
}

fn syndrome_elements(fq: &Fq, alpha: &FqElem, word_bits: &[u64], count: usize) -> PyResult<Vec<FqElem>> {
    let mut out = Vec::with_capacity(count);
    for i in 1..=count {
        let root = fq.pow(alpha, i as i128)?;
        out.push(eval_binary_poly_at(fq, word_bits, &root)?);
    }
    Ok(out)
}

fn berlekamp_massey(fq: &Fq, syndromes: &[FqElem]) -> PyResult<Vec<FqElem>> {
    let zero = fq.zero();
    let one = fq.one();

    let mut c = vec![one.clone()];
    let mut b = vec![one];
    let mut l: usize = 0;
    let mut m: usize = 1;
    let mut bb = fq.one();

    for n in 0..syndromes.len() {
        let mut d = syndromes[n].clone();
        for i in 1..=l {
            let ci = fq_poly_get(&c, i, &zero);
            let si = &syndromes[n - i];
            let term = fq.mul(ci, si)?;
            d = fq.add(&d, &term)?;
        }

        if d.is_zero() {
            m += 1;
            continue;
        }

        let factor = fq.mul(&d, &fq.inv(&bb)?)?;
        let correction = fq_poly_scale_shift(fq, &b, &factor, m)?;
        let t = fq_poly_add(fq, &c, &correction)?;

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

fn chien_search(fq: &Fq, alpha: &FqElem, n: u64, sigma: &[FqElem]) -> PyResult<Vec<u64>> {
    let mut positions = Vec::new();
    for j in 0..n {
        let x = fq.pow(alpha, -(j as i128))?;
        if fq_poly_eval(fq, sigma, &x)?.is_zero() {
            positions.push(j);
        }
    }
    Ok(positions)
}

fn fq_poly_eval(fq: &Fq, coeffs: &[FqElem], x: &FqElem) -> PyResult<FqElem> {
    let mut acc = fq.zero();
    for coeff in coeffs.iter().rev() {
        acc = fq.mul(&acc, x)?;
        acc = fq.add(&acc, coeff)?;
    }
    Ok(acc)
}

fn fq_poly_add(fq: &Fq, a: &[FqElem], b: &[FqElem]) -> PyResult<Vec<FqElem>> {
    let n = a.len().max(b.len());
    let zero = fq.zero();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let ai = fq_poly_get(a, i, &zero);
        let bi = fq_poly_get(b, i, &zero);
        out.push(fq.add(ai, bi)?);
    }
    Ok(fq_poly_trim(out))
}

fn fq_poly_scale_shift(fq: &Fq, poly: &[FqElem], scale: &FqElem, shift: usize) -> PyResult<Vec<FqElem>> {
    let mut out = vec![fq.zero(); poly.len() + shift];
    for (i, coeff) in poly.iter().enumerate() {
        out[i + shift] = fq.mul(coeff, scale)?;
    }
    Ok(fq_poly_trim(out))
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

fn fq_poly_get<'a>(coeffs: &'a [FqElem], idx: usize, zero: &'a FqElem) -> &'a FqElem {
    coeffs.get(idx).unwrap_or(zero)
}

fn binary_poly_from_bits(bits: &[u64]) -> PyResult<PolyFp> {
    if bits.iter().any(|&b| b > 1) {
        return Err(PyValueError::new_err(
            "binary words must contain only 0/1 coefficients",
        ));
    }
    PolyFp::new(2, bits.iter().map(|&b| b as i128).collect())
}

fn bits_from_binary_poly(poly: &PolyFp) -> PyResult<Vec<u64>> {
    poly.coeffs()
        .into_iter()
        .map(|c| match c {
            0 => Ok(0),
            1 => Ok(1),
            _ => Err(PyValueError::new_err(
                "polynomial over F2 has non-binary coefficient",
            )),
        })
        .collect()
}

fn pad_bits(mut bits: Vec<u64>, len: usize) -> Vec<u64> {
    bits.resize(len, 0);
    bits
}

fn shift_poly(poly: &PolyFp, shift: usize) -> PyResult<PolyFp> {
    if poly.is_zero() {
        return PolyFp::new(poly.p(), vec![]);
    }
    let mut coeffs = vec![0i128; shift];
    coeffs.extend(poly.coeffs().into_iter().map(|c| c as i128));
    PolyFp::new(poly.p(), coeffs)
}

fn bits_from_mask(mask: usize, len: usize) -> Vec<u64> {
    (0..len).map(|i| ((mask >> i) & 1) as u64).collect()
}

fn hamming_weight(bits: &[u64]) -> usize {
    bits.iter().filter(|&&b| b != 0).count()
}

#[cfg(test)]
fn multiply_row_vector_matrix(vec: &[u64], matrix: &[Vec<u64>]) -> Vec<u64> {
    if matrix.is_empty() {
        return vec![];
    }
    let cols = matrix[0].len();
    let mut out = vec![0; cols];
    for (i, &bit) in vec.iter().enumerate() {
        if bit == 0 {
            continue;
        }
        for j in 0..cols {
            out[j] ^= matrix[i][j];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hamming_7_4_generator_is_expected() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        assert_eq!(code.length(), 7);
        assert_eq!(code.dimension(), 4);
        assert_eq!(code.generator_poly().coeffs(), vec![1, 1, 0, 1]);
    }

    #[test]
    fn bch_encode_produces_divisible_codeword() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let msg = vec![1, 0, 1, 1];
        let word = code.encode(msg.clone()).unwrap();
        assert_eq!(word.len(), 7);
        assert!(code.is_codeword(word.clone()).unwrap());
        assert!(code.parity_check(word).unwrap());
        assert_eq!(code.extract_message(code.encode(msg.clone()).unwrap()).unwrap(), msg);
    }

    #[test]
    fn systematic_encoding_round_trips_message() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let msg = vec![1, 0, 1, 1];
        let word = code.encode_systematic(msg.clone()).unwrap();
        assert!(code.is_codeword(word.clone()).unwrap());
        assert_eq!(code.extract_systematic_message(word).unwrap(), msg);
    }

    #[test]
    fn shortened_systematic_round_trips_message() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        assert_eq!(code.shortened_parameters(3).unwrap(), (12, 4));
        let msg = vec![1, 0, 1, 1];
        let word = code.encode_shortened_systematic(msg.clone(), 3).unwrap();
        assert_eq!(word.len(), 12);
        assert_eq!(
            code.extract_shortened_systematic_message(word, 3).unwrap(),
            msg
        );
    }

    #[test]
    fn parity_check_detects_single_bit_error() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let mut word = code.encode(vec![1, 1, 0, 1]).unwrap();
        word[2] ^= 1;
        assert!(!code.is_codeword(word.clone()).unwrap());
        assert!(!code.parity_check(word.clone()).unwrap());
        assert!(code.syndromes(word).unwrap().iter().any(|s| !s.is_empty()));
    }

    #[test]
    fn binary_bch_15_7_has_known_generator_degree() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        assert_eq!(code.length(), 15);
        assert_eq!(code.generator_degree(), 8);
        assert_eq!(code.dimension(), 7);
        assert_eq!(code.generator_poly().coeffs().len(), 9);
    }

    #[test]
    fn hamming_decoder_corrects_single_error() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let msg = vec![1, 0, 1, 1];
        let sent = code.encode(msg.clone()).unwrap();
        let mut received = sent.clone();
        received[4] ^= 1;

        let corrected = code.decode(received).unwrap();
        assert_eq!(corrected, sent);
        assert!(code.parity_check(corrected).unwrap());
        assert_eq!(code.decode_message(vec![1, 1, 1, 1, 0, 1, 1]).unwrap(), msg);
    }

    #[test]
    fn systematic_decoder_recovers_message_after_error_correction() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let msg = vec![1, 0, 1, 1];
        let mut received = code.encode_systematic(msg.clone()).unwrap();
        received[5] ^= 1;
        assert_eq!(code.decode_systematic_message(received).unwrap(), msg);
    }

    #[test]
    fn shortened_decoder_recovers_message_after_error_correction() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        let msg = vec![1, 0, 1, 1];
        let mut received = code.encode_shortened_systematic(msg.clone(), 3).unwrap();
        received[2] ^= 1;
        received[8] ^= 1;
        assert_eq!(
            code.decode_shortened_systematic_message(received, 3).unwrap(),
            msg
        );
    }

    #[test]
    fn hamming_decoder_corrects_all_single_bit_errors_exhaustively() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        for mask in 0usize..(1 << code.dimension()) {
            let msg = bits_from_mask(mask, code.dimension());
            let codeword = code.encode_systematic(msg.clone()).unwrap();
            for pos in 0..(code.length() as usize) {
                let mut received = codeword.clone();
                received[pos] ^= 1;
                assert_eq!(code.decode(received.clone()).unwrap(), codeword);
                assert_eq!(code.decode_systematic_message(received).unwrap(), msg);
            }
        }
    }

    #[test]
    fn shortened_bch_decoding_matches_full_length_padding() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        let shorten_by = 3usize;
        let k_short = code.shortened_parameters(shorten_by).unwrap().1;
        for mask in 0usize..(1 << k_short) {
            let msg = bits_from_mask(mask, k_short);
            let short_codeword = code.encode_shortened_systematic(msg.clone(), shorten_by).unwrap();

            let mut short_received = short_codeword.clone();
            short_received[0] ^= 1;
            short_received[5] ^= 1;

            let short_decoded = code
                .decode_shortened_systematic_message(short_received.clone(), shorten_by)
                .unwrap();

            let mut padded_received = short_received;
            padded_received.resize(code.length() as usize, 0);
            let full_decoded = code.decode_systematic_message(padded_received).unwrap();
            assert_eq!(short_decoded, msg);
            assert_eq!(full_decoded[..k_short], msg);
        }
    }

    #[test]
    fn bch_15_7_5_decoder_corrects_two_errors() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        let msg = vec![1, 0, 1, 1, 0, 1, 1];
        let sent = code.encode(msg.clone()).unwrap();
        let mut received = sent.clone();
        received[2] ^= 1;
        received[11] ^= 1;

        let corrected = code.decode(received).unwrap();
        assert_eq!(corrected, sent);
        assert_eq!(code.extract_message(sent).unwrap(), msg);
    }

    #[test]
    fn extract_message_rejects_non_codeword() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        assert!(code.extract_message(vec![1, 0, 0, 0, 0, 0, 0]).is_err());
    }

    #[test]
    fn decode_returns_codeword_on_three_errors_without_claiming_original() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        let sent = code.encode(vec![1, 0, 1, 1, 0, 1, 1]).unwrap();
        let mut received = sent.clone();
        received[1] ^= 1;
        received[4] ^= 1;
        received[7] ^= 1;

        let decoded = code.decode(received).unwrap();
        assert!(code.is_codeword(decoded).unwrap());
    }

    #[test]
    fn generator_and_parity_check_matrices_are_orthogonal() {
        let code = BinaryBchCode::new(4, vec![1, 1, 0, 0, 1], 5).unwrap();
        let g = code.generator_matrix().unwrap();
        let h = code.parity_check_matrix().unwrap();
        let ht = transpose(&h);
        for row in &g {
            let syndrome = multiply_row_vector_matrix(row, &ht);
            assert!(syndrome.iter().all(|&b| b == 0));
        }
    }

    #[test]
    fn hamming_code_enumeration_and_distance_are_correct() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let codewords = code.codewords(Some(8)).unwrap();
        assert_eq!(codewords.len(), 16);
        assert_eq!(code.minimum_distance(Some(8)).unwrap(), 3);
    }

    #[test]
    fn hamming_weight_distribution_matches_known_values() {
        let code = BinaryBchCode::new(3, vec![1, 1, 0, 1], 3).unwrap();
        let dist = code.weight_distribution(Some(8)).unwrap();
        assert_eq!(dist[0], 1);
        assert_eq!(dist[3], 7);
        assert_eq!(dist[4], 7);
        assert_eq!(dist[7], 1);
        assert_eq!(dist.iter().sum::<u64>(), 16);
    }

    fn transpose(matrix: &[Vec<u64>]) -> Vec<Vec<u64>> {
        if matrix.is_empty() {
            return vec![];
        }
        let rows = matrix.len();
        let cols = matrix[0].len();
        let mut out = vec![vec![0; rows]; cols];
        for (i, row) in matrix.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                out[j][i] = val;
            }
        }
        out
    }
}

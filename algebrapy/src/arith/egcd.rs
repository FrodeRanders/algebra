/// Extended Euclidean algorithm.
/// Returns (g, x, y) such that a*x + b*y = g = gcd(a,b)
pub fn egcd_i128(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a.abs(), a.signum(), 0)
    } else {
        let (g, x1, y1) = egcd_i128(b, a % b);
        let x = y1;
        let y = x1 - (a / b) * y1;
        (g, x, y)
    }
}

/// Modular inverse of a mod m (for gcd(a,m)=1).
/// Returns value in [0..m-1] as i128.
pub fn inv_mod_i128(a: i128, m: i128) -> Option<i128> {
    let (g, x, _) = egcd_i128(a, m);
    if g != 1 {
        return None;
    }
    let mut r = x % m;
    if r < 0 {
        r += m;
    }
    Some(r)
}

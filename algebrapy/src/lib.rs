use pyo3::prelude::*;

pub mod arith;
pub mod field;
pub mod group;

use field::fp::{Fp, FpElem};
use field::poly_fp::PolyFp;
use field::fq::{Fq, FqElem};
use group::perm::{Perm, Sn};

/// A Python module implemented in Rust.
#[pymodule]
mod algebrapy {
    use super::*;

    fn register_classes(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Fp>()?;
        m.add_class::<FpElem>()?;
        m.add_class::<PolyFp>()?;
        m.add_class::<Fq>()?;
        m.add_class::<FqElem>()?;
        m.add_class::<Perm>()?;
        m.add_class::<Sn>()?;
        Ok(())
    }

    #[pymodule_init]
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        register_classes(m)?;
        Ok(())
    }
}

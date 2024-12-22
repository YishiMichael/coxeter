use alga::general as alg;

use super::graded_algebra::GradedAlgebra;
use super::graded_algebra::GradedAlgebraMeta;

#[derive(Clone)]
pub struct PolynomialMeta<R> {
    phantom: std::marker::PhantomData<R>,
}

impl<R> GradedAlgebraMeta for PolynomialMeta<R>
where
    R: alg::RingCommutative,
{
    type Ring = R;
    type Grade = u32;
}

pub type Polynomial<R> = GradedAlgebra<PolynomialMeta<R>>;

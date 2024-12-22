use alga::general as alg;

use super::graded_algebra::GradedAlgebra;
use super::graded_algebra::GradedAlgebraMeta;

#[derive(Clone, Debug)]
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

impl<R: alg::RingCommutative> Polynomial<R> {
    pub fn monomial(coeff: R, exponent: u32) -> Self {
        Self::embed(coeff, exponent)
    }

    pub fn eval(self, value: R) -> R {
        self.into_terms()
            .map(|(coeff, exponent)| {
                std::iter::repeat(value.clone())
                    .take(exponent as usize)
                    .fold(coeff, R::mul)
            })
            .fold(R::zero(), R::add)
    }
}

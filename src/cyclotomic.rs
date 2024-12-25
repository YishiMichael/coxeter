use alga::general as alg;
use itertools::Itertools;

use super::graded_algebra::GradedAlgebra;
use super::graded_algebra::GradedAlgebraMeta;
use super::prufer_monoid::PruferMonoid;

#[derive(Clone, Debug)]
pub struct CyclotomicMeta<R> {
    phantom: std::marker::PhantomData<R>,
}

impl<R> GradedAlgebraMeta for CyclotomicMeta<R>
where
    R: alg::RingCommutative,
{
    type Ring = R;
    type Grade = PruferMonoid<u64>;

    fn annihilate(element: GradedAlgebra<Self>) -> GradedAlgebra<Self> {
        fn prime_factors(n: u64) -> impl Iterator<Item = u64> {
            std::iter::repeat(()).scan((n, 2), |(n, p), _| {
                if *n == 1 {
                    None
                } else {
                    *p = (*p..).find(|p| *n % p == 0).unwrap();
                    *n /= *p;
                    Some(*p)
                }
            })
        }

        prime_factors(
            element
                .clone()
                .into_terms()
                .map(|(_, grade)| grade.denom().clone())
                .fold(1, |common_multiple, denominator| {
                    num::Integer::lcm(&common_multiple, &denominator)
                }),
        )
        .dedup()
        .map(|factor| {
            (1..factor)
                .map(|i| {
                    Cyclotomic::root_of_unity(factor, 0) - Cyclotomic::root_of_unity(factor, i)
                })
                .fold(num::Zero::zero(), std::ops::Add::add)
        })
        .fold(element, std::ops::Mul::mul)
    }
}

pub type Cyclotomic<R> = GradedAlgebra<CyclotomicMeta<R>>;

impl<R: alg::RingCommutative> Cyclotomic<R> {
    pub fn root_of_unity(order: u64, exponent: u64) -> Self {
        Self::embed(
            R::one(),
            PruferMonoid::new(num_rational::Ratio::new(exponent, order)),
        )
    }
}

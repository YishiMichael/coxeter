use alga::general as alg;
use itertools::Itertools;
use num::One;

use super::graded_algebra::GradedAlgebra;
use super::graded_algebra::GradedAlgebraMeta;
use super::prufer_group::PruferGroup;

#[derive(Clone, Debug)]
pub struct CyclotomicMeta<R> {
    phantom: std::marker::PhantomData<R>,
}

impl<R> GradedAlgebraMeta for CyclotomicMeta<R>
where
    R: alg::RingCommutative,
{
    type Ring = R;
    type Grade = PruferGroup<u32>;

    fn annihilate(element: GradedAlgebra<Self>) -> GradedAlgebra<Self> {
        fn prime_factors(n: u32) -> impl Iterator<Item = u32> {
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
                .map(|(_, grade)| grade.denominator)
                .fold(1, |common_multiple, denominator| {
                    num::Integer::lcm(&common_multiple, &denominator)
                }),
        )
        .dedup()
        .map(|factor| {
            (1..factor)
                .map(|i| Cyclotomic::one() - Cyclotomic::root_of_unity(factor, i))
                .fold(
                    <Cyclotomic<R> as num::Zero>::zero(),
                    <Cyclotomic<R> as std::ops::Add>::add,
                )
        })
        .fold(element, <Cyclotomic<R> as std::ops::Mul>::mul)
    }
}

pub type Cyclotomic<R> = GradedAlgebra<CyclotomicMeta<R>>;

impl<R: alg::RingCommutative> Cyclotomic<R> {
    pub fn root_of_unity(order: u32, exponent: u32) -> Self {
        Self::embed(R::one(), PruferGroup::new(order, exponent))
    }
}

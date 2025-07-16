use itertools::Itertools;

use super::alg::*;
use super::graded_algebra::GradedAlgebra;
use super::prufer_group::PruferGroup;

#[derive(Clone)]
pub struct Cyclotomic<C, U>(GradedAlgebra<C, PruferGroup<U>>);

impl<C, U> Cyclotomic<C, U>
where
    C: Ring,
    U: num::Integer + num::Unsigned,
{
    pub fn root_of_unity(order: U, exponent: U) -> Self {
        Self(GradedAlgebra::from_terms(std::iter::once((
            C::one(),
            PruferGroup::new(exponent, order),
        ))))
    }
}

impl<C, U> PartialEq for Cyclotomic<C, U>
where
    C: Ring + Clone,
    U: num::Integer + num::Unsigned + Clone,
{
    fn eq(&self, rhs: &Self) -> bool {
        self.clone().sub(rhs.clone()).is_zero()
    }
}

impl<C, U> AdditiveMagma for Cyclotomic<C, U>
where
    C: Ring + Clone,
    U: num::Integer + num::Unsigned + Clone,
{
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(GradedAlgebra::sum(iter.map(|element| element.0)))
    }
}

impl<C, U> AdditiveIdentity for Cyclotomic<C, U>
where
    C: Ring + Clone,
    U: num::Integer + num::Unsigned + Clone,
{
    #[inline]
    fn zero() -> Self {
        Self(GradedAlgebra::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        let order = self
            .0
            .get_terms()
            .map(|(_, exponent)| exponent.denom())
            .fold(U::one(), |denom_0, denom_1| denom_0.lcm(denom_1));
        let balancer = GradedAlgebra::prod(
            std::iter::successors(Some(U::one() + U::one()), |p| Some(p.clone().add(U::one())))
                .scan(order, |n, mut p| {
                    (!n.is_one()).then(|| {
                        while !n.is_multiple_of(&p) {
                            p.inc();
                        }
                        *n = n.div_floor(&p);
                        p
                    })
                })
                .dedup()
                .map(|p| {
                    GradedAlgebra::sum(
                        std::iter::successors(Some(U::zero()), |i| Some(i.clone().add(U::one())))
                            .take_while(|i| *i < p)
                            .map(|i| {
                                Self::root_of_unity(p.clone(), U::zero())
                                    .0
                                    .sub(Self::root_of_unity(p.clone(), i).0)
                            }),
                    )
                }),
        );
        self.0.clone().mul(balancer).is_zero()
    }
}

impl<C, U> AdditiveInverse for Cyclotomic<C, U>
where
    C: Ring + Clone,
    U: num::Integer + num::Unsigned + Clone,
{
    #[inline]
    fn neg(self) -> Self {
        Self(self.0.neg())
    }
}

impl<C, U> MultiplicativeMagma for Cyclotomic<C, U>
where
    C: Ring + Clone,
    U: num::Integer + num::Unsigned + Clone,
{
    #[inline]
    fn prod<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(GradedAlgebra::prod(iter.map(|element| element.0)))
    }
}

impl<C, U> MultiplicativeIdentity for Cyclotomic<C, U>
where
    C: Ring + Clone,
    U: num::Integer + num::Unsigned + Clone,
{
    #[inline]
    fn one() -> Self {
        Self::root_of_unity(U::one(), U::zero())
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.eq(&Self::one())
    }
}

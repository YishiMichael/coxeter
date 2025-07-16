use itertools::Itertools;

use super::alg::*;
// use super::big_uint::BigUint;
use super::graded_algebra::GradedAlgebra;
use super::prufer_group::PruferGroup;
// use super::rational::Rational;

#[derive(Clone, Debug)]
pub struct Cyclotomic<C, U>(GradedAlgebra<C, PruferGroup<U>>);

impl<C, U> Cyclotomic<C, U>
where
    C: Ring,
    U: num::Integer + num::Unsigned,
    // PruferGroup<U>: AdditiveGroup,
{
    pub fn root_of_unity(order: U, exponent: U) -> Self {
        // let order = order.into();
        // let exponent = exponent.into();
        // Self(GradedAlgebra::mul(
        //     Self::root_of_unity_graded(order.clone(), exponent),
        //     Self::balancer(order.prime_factors()),
        // ))
        Self(GradedAlgebra::from_terms(std::iter::once((
            C::one(),
            PruferGroup::new(exponent, order),
        ))))
    }

    // fn root_of_unity_graded(
    //     order: BigUint,
    //     exponent: BigUint,
    // ) -> GradedAlgebra<Rational, PruferGroup> {
    //     GradedAlgebra::embed(
    //         Rational::one(),
    //         PruferGroup::from((exponent.into(), order.into())),
    //     )
    // }

    // fn balancer<I: Iterator<Item = BigUint>>(
    //     prime_factors: I,
    // ) -> GradedAlgebra<BigRational, PruferGroup> {
    //     GradedAlgebra::product(prime_factors.map(|factor| {
    //         GradedAlgebra::mul(
    //             GradedAlgebra::embed_scalar((BigUint::one(), factor.clone()).into()),
    //             GradedAlgebra::sum(factor.clone().range().map(|i| {
    //                 GradedAlgebra::sub(
    //                     Self::root_of_unity_graded(factor.clone(), BigUint::zero()),
    //                     Self::root_of_unity_graded(factor.clone(), i),
    //                 )
    //             })),
    //         )
    //     }))
    // }

    // fn raise_order(
    //     value: GradedAlgebra<BigRational, PruferGroup>,
    //     low_order: BigUint,
    //     high_order: BigUint,
    // ) -> GradedAlgebra<BigRational, PruferGroup> {
    //     if low_order == high_order {
    //         value
    //     } else {
    //         println!("##1, {value:?} | {low_order:?} -> {high_order:?}");
    //         let prime_factors = low_order.prime_factors().collect::<Vec<_>>();
    //         println!("##2");
    //         let additional_prime_factors = high_order
    //             .clone()
    //             .prime_factors()
    //             .dedup()
    //             .filter(|order_factor| prime_factors.iter().all(|factor| factor != order_factor))
    //             .collect::<Vec<_>>();
    //         println!("##3");
    //         if additional_prime_factors.is_empty() {
    //             value
    //         } else {
    //             println!("##4, {additional_prime_factors:?}");
    //             GradedAlgebra::mul(value, Self::balancer(additional_prime_factors.into_iter()))
    //         }
    //     }
    // }
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
        // let elements = iter.collect::<Vec<_>>();
        // let order = BigUint::lcm(elements.iter().map(|element| element.order.clone()));
        // Self {
        //     order: order.clone(),
        //     value: GradedAlgebra::product(
        //         elements
        //             .into_iter()
        //             .map(|element| Self::raise_order(element.value, element.order, order.clone())),
        //     ),
        // }
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

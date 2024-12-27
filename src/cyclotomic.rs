use itertools::Itertools;

use super::alg;
use super::big_rational::BigRational;
use super::big_uint::BigUint;
use super::graded_algebra::GradedAlgebra;
use super::prufer_group::PruferGroup;

#[derive(Clone, Debug)]
pub struct Cyclotomic {
    order: BigUint,
    value: GradedAlgebra<BigRational, PruferGroup>,
}

impl Cyclotomic {
    pub fn root_of_unity<T>(order: T, exponent: T) -> Self
    where
        T: Into<BigUint>,
    {
        let order = order.into();
        let exponent = exponent.into();
        Self {
            order: order.clone(),
            value: alg::MultiplicativeMagma::mul(
                Self::root_of_unity_raw_value(order.clone(), exponent),
                Self::balancer(order.prime_factors()),
            ),
        }
    }

    fn root_of_unity_raw_value(
        order: BigUint,
        exponent: BigUint,
    ) -> GradedAlgebra<BigRational, PruferGroup> {
        GradedAlgebra::embed(
            alg::MultiplicativeIdentity::one(),
            PruferGroup::from((exponent.into(), order.into())),
        )
    }

    fn balancer<I: Iterator<Item = BigUint>>(
        prime_factors: I,
    ) -> GradedAlgebra<BigRational, PruferGroup> {
        alg::MultiplicativeMagma::product(prime_factors.map(|factor| {
            alg::MultiplicativeMagma::mul(
                GradedAlgebra::embed_scalar(
                    (alg::MultiplicativeIdentity::one(), factor.clone()).into(),
                ),
                alg::AdditiveMagma::sum(factor.clone().range().map(|i| {
                    alg::AdditiveMagma::sub(
                        Self::root_of_unity_raw_value(
                            factor.clone(),
                            alg::AdditiveIdentity::zero(),
                        ),
                        Self::root_of_unity_raw_value(factor.clone(), i),
                    )
                })),
            )
        }))
    }

    fn raise_order(
        value: GradedAlgebra<BigRational, PruferGroup>,
        low_order: BigUint,
        high_order: BigUint,
    ) -> GradedAlgebra<BigRational, PruferGroup> {
        if low_order == high_order {
            value
        } else {
            let prime_factors = low_order.prime_factors().collect::<Vec<_>>();
            let additional_prime_factors = high_order
                .clone()
                .prime_factors()
                .dedup()
                .filter(|order_factor| prime_factors.iter().all(|factor| factor != order_factor))
                .collect::<Vec<_>>();
            if additional_prime_factors.is_empty() {
                value
            } else {
                alg::MultiplicativeMagma::mul(
                    value,
                    Self::balancer(additional_prime_factors.into_iter()),
                )
            }
        }
    }
}

impl PartialEq for Cyclotomic {
    fn eq(&self, rhs: &Self) -> bool {
        alg::AdditiveIdentity::is_zero(&alg::AdditiveMagma::sub(self.clone(), rhs.clone()))
    }
}

impl alg::AdditiveMagma for Cyclotomic {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let elements = iter.collect::<Vec<_>>();
        let order = BigUint::lcm(elements.iter().map(|element| element.order.clone()));
        Self {
            order: order.clone(),
            value: alg::AdditiveMagma::sum(
                elements
                    .into_iter()
                    .map(|element| Self::raise_order(element.value, element.order, order.clone())),
            ),
        }
    }
}

impl alg::AdditiveIdentity for Cyclotomic {
    #[inline]
    fn zero() -> Self {
        Self {
            order: alg::MultiplicativeIdentity::one(),
            value: alg::AdditiveIdentity::zero(),
        }
    }

    #[inline]
    fn is_zero(&self) -> bool {
        alg::AdditiveIdentity::is_zero(&self.value)
    }
}

impl alg::AdditiveInverse for Cyclotomic {
    #[inline]
    fn neg(self) -> Self {
        Self {
            order: self.order,
            value: self.value.neg(),
        }
    }
}

impl alg::MultiplicativeMagma for Cyclotomic {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let elements = iter.collect::<Vec<_>>();
        let order = BigUint::lcm(elements.iter().map(|element| element.order.clone()));
        Self {
            order: order.clone(),
            value: alg::MultiplicativeMagma::product(
                elements
                    .into_iter()
                    .map(|element| Self::raise_order(element.value, element.order, order.clone())),
            ),
        }
    }
}

impl alg::MultiplicativeIdentity for Cyclotomic {
    #[inline]
    fn one() -> Self {
        Self::root_of_unity(1, 0)
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.eq(&Self::one())
    }
}

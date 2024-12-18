use std::collections::HashMap;

use abstalg as alg;
use itertools::Itertools;

lazy_static::lazy_static! {
    static ref POLYNOMIAL_ALGEBRA: alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>
        = alg::PolynomialAlgebra::new(alg::ReducedFractions::new(alg::I32));

    static ref CYCLOTOMIC_FIELDS_CACHE: HashMap<usize, alg::QuotientField<alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>>>
        = HashMap::new();
}

fn cyclotomic_field(
    order: usize,
) -> &'static alg::QuotientField<alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>>
{
    CYCLOTOMIC_FIELDS_CACHE.entry(order).or_insert_with(|| {
        // Construct the cyclotomic polynomial with Möbius inversion formula.
        // See https://en.wikipedia.org/wiki/Cyclotomic_polynomial
        let modulo = cyclotomic_polynomial(POLYNOMIAL_ALGEBRA, order);
        alg::QuotientField::new(POLYNOMIAL_ALGEBRA.clone(), modulo)
    })
}

fn cyclotomic_polynomial(
    polynomial_algebra: &alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>,
    order: usize,
) -> <alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>> as alg::Domain>::Elem {
    let (numerator, denominator) = prime_factors(order).dedup().powerset().fold(
        (
            alg::Monoid::one(polynomial_algebra),
            alg::Monoid::one(polynomial_algebra),
        ),
        |(numerator, denominator), factors| {
            let is_denominator = factors.len() % 2 == 1;
            let degree = order / factors.into_iter().product::<usize>();
            // x^d - 1
            let polynomial_factor = alg::AbelianGroup::sub(
                polynomial_algebra,
                &alg::CommuntativeMonoid::times(
                    polynomial_algebra,
                    degree,
                    &vec![
                        alg::CommuntativeMonoid::zero(polynomial_algebra.base()),
                        alg::Monoid::one(polynomial_algebra.base()),
                    ],
                ),
                &alg::Monoid::one(polynomial_algebra),
            );
            if is_denominator {
                (
                    numerator,
                    alg::Semigroup::mul(polynomial_algebra, &polynomial_factor, &denominator),
                )
            } else {
                (
                    alg::Semigroup::mul(polynomial_algebra, &polynomial_factor, &numerator),
                    denominator,
                )
            }
        },
    );
    alg::EuclideanDomain::quo(polynomial_algebra, &numerator, &denominator)
}

fn prime_factors(n: usize) -> impl Iterator<Item = usize> {
    std::iter::repeat(()).scan((n, 2usize), |(n, p), _| {
        if *n == 1 {
            None
        } else {
            *p = (*p..).find(|p| *n % p == 0).unwrap();
            *n /= *p;
            Some(*p)
        }
    })
}

// #[derive(Clone, Debug)]
// pub struct CyclotomicFactory<A>
// where
//     A: AbelianGroup,
// {
//     base: A,
//     len: usize,
// }

// impl<A> Cyclotomic<A>
// where
//     A: AbelianGroup,
// {
//     pub fn new(base: A, len: usize) -> Self {
//         Self { base, len }
//     }
// }

pub trait GroupElement:
    Clone
    + std::ops::Add<Output = Self>
    + std::ops::AddAssign
    + std::ops::Sub<Output = Self>
    + std::ops::SubAssign
    + std::ops::Neg<Output = Self>
    + num::Zero
{
}

pub trait RingElement:
    GroupElement + std::ops::Mul<Output = Self> + std::ops::MulAssign + num::One
{
}

impl GroupElement for i32 {}
impl RingElement for i32 {}

// Ring structure R[zeta_n]
// Basis of coefficients: see Thomas, Integral Bases for Subfields of Cyclotomic Fields
#[derive(Clone)]
pub struct CyclotomicNumber<C>(Vec<C>);

impl<C: GroupElement> CyclotomicNumber<C> {
    pub fn order(&self) -> usize {
        self.0.len()
    }

    fn zipped(lhs: Vec<C>, rhs: Vec<C>) -> Vec<(C, C)> {
        let lhs_order = lhs.len();
        let rhs_order = rhs.len();
        let order = num::Integer::lcm(&lhs_order, &rhs_order);
        let lhs_scalar = order / lhs_order;
        let rhs_scalar = order / rhs_order;
        (0..order)
            .map(|index| {
                (
                    if index % lhs_scalar == 0 {
                        lhs[index / lhs_scalar].clone()
                    } else {
                        C::zero()
                    },
                    if index % rhs_scalar == 0 {
                        rhs[index / rhs_scalar].clone()
                    } else {
                        C::zero()
                    },
                )
            })
            .collect()
    }
}

impl<C: RingElement> CyclotomicNumber<C> {
    pub fn root_of_unity(order: usize, exponent: usize) -> Self {
        let mut number = Self(vec![C::zero(); order]);
        number.0[exponent] = C::one();
        number
    }
}

impl<C: GroupElement> std::ops::Add for CyclotomicNumber<C> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<C: GroupElement> std::ops::AddAssign for CyclotomicNumber<C> {
    fn add_assign(&mut self, rhs: Self) {
        let zipped = Self::zipped(std::mem::replace(&mut self.0, Vec::new()), rhs.0);
        self.0
            .extend(zipped.into_iter().map(|(lhs, rhs)| lhs + rhs));
    }
}

impl<C: GroupElement> std::ops::Sub for CyclotomicNumber<C> {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<C: GroupElement> std::ops::SubAssign for CyclotomicNumber<C> {
    fn sub_assign(&mut self, rhs: Self) {
        let zipped = Self::zipped(std::mem::replace(&mut self.0, Vec::new()), rhs.0);
        self.0
            .extend(zipped.into_iter().map(|(lhs, rhs)| lhs - rhs));
    }
}

impl<C: GroupElement> std::ops::Neg for CyclotomicNumber<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.into_iter().map(std::ops::Neg::neg).collect())
    }
}

impl<C: GroupElement> num::Zero for CyclotomicNumber<C> {
    fn zero() -> Self {
        Self(vec![C::zero(); 1])
    }

    fn is_zero(&self) -> bool {
        todo!()
    }
}

impl<C: RingElement> std::ops::Mul for CyclotomicNumber<C> {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<C: RingElement> std::ops::MulAssign for CyclotomicNumber<C> {
    fn mul_assign(&mut self, rhs: Self) {
        fn cross_inner_product<C: RingElement>(zipped: &[(C, C)]) -> C {
            zipped
                .iter()
                .map(|(elem_0, _)| elem_0)
                .zip(zipped.iter().map(|(_, elem_1)| elem_1).rev())
                .map(|(elem_0, elem_1)| elem_0.clone() * elem_1.clone())
                .fold(C::zero(), C::add)
        }
        let zipped = Self::zipped(std::mem::replace(&mut self.0, Vec::new()), rhs.0);
        self.0.extend((1..=zipped.len()).map(|index| {
            cross_inner_product(&zipped[..index]) + cross_inner_product(&zipped[index..])
        }))
    }
}

impl<C: RingElement> num::One for CyclotomicNumber<C> {
    fn one() -> Self {
        Self::root_of_unity(1, 0)
    }
}

impl<C: GroupElement> GroupElement for CyclotomicNumber<C> {}
impl<C: RingElement> RingElement for CyclotomicNumber<C> {}

use alga::general as alg;

// use std::collections::HashMap;

// lazy_static::lazy_static! {
//     static ref POLYNOMIAL_ALGEBRA: alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>
//         = alg::PolynomialAlgebra::new(alg::ReducedFractions::new(alg::I32));

//     static ref CYCLOTOMIC_FIELDS_CACHE: HashMap<usize, alg::QuotientField<alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>>>
//         = HashMap::new();
// }

// fn cyclotomic_field(
//     order: usize,
// ) -> &'static alg::QuotientField<alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>>
// {
//     CYCLOTOMIC_FIELDS_CACHE.entry(order).or_insert_with(|| {
//         // Construct the cyclotomic polynomial with Möbius inversion formula.
//         // See https://en.wikipedia.org/wiki/Cyclotomic_polynomial
//         let modulo = cyclotomic_polynomial(POLYNOMIAL_ALGEBRA, order);
//         alg::QuotientField::new(POLYNOMIAL_ALGEBRA.clone(), modulo)
//     })
// }

// fn cyclotomic_polynomial(
//     polynomial_algebra: &alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>>,
//     order: usize,
// ) -> <alg::PolynomialAlgebra<alg::ReducedFractions<alg::CheckedInts<i32>>> as alg::Domain>::Elem {
//     let (numerator, denominator) = prime_factors(order).dedup().powerset().fold(
//         (
//             alg::Monoid::one(polynomial_algebra),
//             alg::Monoid::one(polynomial_algebra),
//         ),
//         |(numerator, denominator), factors| {
//             let is_denominator = factors.len() % 2 == 1;
//             let degree = order / factors.into_iter().product::<usize>();
//             // x^d - 1
//             let polynomial_factor = alg::AbelianGroup::sub(
//                 polynomial_algebra,
//                 &alg::CommuntativeMonoid::times(
//                     polynomial_algebra,
//                     degree,
//                     &vec![
//                         alg::CommuntativeMonoid::zero(polynomial_algebra.base()),
//                         alg::Monoid::one(polynomial_algebra.base()),
//                     ],
//                 ),
//                 &alg::Monoid::one(polynomial_algebra),
//             );
//             if is_denominator {
//                 (
//                     numerator,
//                     alg::Semigroup::mul(polynomial_algebra, &polynomial_factor, &denominator),
//                 )
//             } else {
//                 (
//                     alg::Semigroup::mul(polynomial_algebra, &polynomial_factor, &numerator),
//                     denominator,
//                 )
//             }
//         },
//     );
//     alg::EuclideanDomain::quo(polynomial_algebra, &numerator, &denominator)
// }

// fn prime_factors(n: usize) -> impl Iterator<Item = usize> {
//     std::iter::repeat(()).scan((n, 2usize), |(n, p), _| {
//         if *n == 1 {
//             None
//         } else {
//             *p = (*p..).find(|p| *n % p == 0).unwrap();
//             *n /= *p;
//             Some(*p)
//         }
//     })
// }

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

// pub trait GroupElement:
//     Clone
//     + std::ops::Add<Output = Self>
//     + std::ops::AddAssign
//     + std::ops::Sub<Output = Self>
//     + std::ops::SubAssign
//     + std::ops::Neg<Output = Self>
//     + num::Zero
// {
// }

// pub trait RingElement:
//     GroupElement + std::ops::Mul<Output = Self> + std::ops::MulAssign + num::One
// {
// }

// impl GroupElement for i32 {}
// impl RingElement for i32 {}

// A ring adjoint with roots of unity.
// The length of vector represents the order.
// Basis of coefficients: see Thomas, Integral Bases for Subfields of Cyclotomic Fields
#[derive(Clone, Debug)]
pub struct Cyclotomic<C>(Vec<C>);

impl<C> Cyclotomic<C> {
    // pub fn order(&self) -> usize {
    //     self.0.len()
    // }

    pub fn root_of_unity(order: usize, exponent: usize) -> Self
    where
        C: alg::Identity<alg::Additive> + alg::Identity<alg::Multiplicative>,
    {
        // let mut number = std::iter::repeat_with(|| <C as alg::Identity<alg::Additive>>::identity()).take(order);
        // number.0[exponent] = <C as alg::Identity<alg::Multiplicative>>::identity();
        // number
    }
}

fn zipped<C: alg::RingCommutative>(lhs: Vec<C>, rhs: Vec<C>) -> Vec<(C, C)> {
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

impl<C: alg::RingCommutative> std::ops::Add for Cyclotomic<C> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<C: alg::RingCommutative> std::ops::AddAssign for Cyclotomic<C> {
    fn add_assign(&mut self, rhs: Self) {
        let zipped = Self::zipped(std::mem::replace(&mut self.0, Vec::new()), rhs.0);
        self.0
            .extend(zipped.into_iter().map(|(lhs, rhs)| lhs + rhs));
    }
}

impl<C: alg::RingCommutative> std::ops::Sub for Cyclotomic<C> {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<C: alg::RingCommutative> std::ops::SubAssign for Cyclotomic<C> {
    fn sub_assign(&mut self, rhs: Self) {
        let zipped = Self::zipped(std::mem::replace(&mut self.0, Vec::new()), rhs.0);
        self.0
            .extend(zipped.into_iter().map(|(lhs, rhs)| lhs - rhs));
    }
}

impl<C: alg::RingCommutative> std::ops::Neg for Cyclotomic<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.into_iter().map(std::ops::Neg::neg).collect())
    }
}

impl<C: alg::RingCommutative> num::Zero for Cyclotomic<C> {
    fn zero() -> Self {
        Self(vec![C::zero(); 1])
    }

    fn is_zero(&self) -> bool {
        todo!()
    }
}

impl<C: alg::RingCommutative> std::ops::Mul for Cyclotomic<C> {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<C: alg::RingCommutative> std::ops::MulAssign for Cyclotomic<C> {
    fn mul_assign(&mut self, rhs: Self) {
        fn cross_inner_product<C: alg::RingCommutative>(zipped: &[(C, C)]) -> C {
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

impl<C: alg::RingCommutative> num::One for Cyclotomic<C> {
    fn one() -> Self {
        Self::root_of_unity(1, 0)
    }
}

////////////////////////
// TODO: specify trait bounds
impl<C: alg::RingCommutative> PartialEq for Cyclotomic<C> {
    fn eq(&self, _rhs: &Self) -> bool {
        todo!()
    }
}

impl<C: alg::RingCommutative> alg::AbstractMagma<alg::Additive> for Cyclotomic<C> {
    fn operate(&self, _rhs: &Self) -> Self {
        todo!()
    }
}

impl<C: alg::RingCommutative> alg::TwoSidedInverse<alg::Additive> for Cyclotomic<C> {
    fn two_sided_inverse(&self) -> Self {
        todo!()
    }
}

impl<C: alg::RingCommutative> alg::Identity<alg::Additive> for Cyclotomic<C> {
    fn identity() -> Self {
        todo!()
    }
}

impl<C: alg::RingCommutative> alg::AbstractQuasigroup<alg::Additive> for Cyclotomic<C> {}
impl<C: alg::RingCommutative> alg::AbstractLoop<alg::Additive> for Cyclotomic<C> {}
impl<C: alg::RingCommutative> alg::AbstractSemigroup<alg::Additive> for Cyclotomic<C> {}
impl<C: alg::RingCommutative> alg::AbstractMonoid<alg::Additive> for Cyclotomic<C> {}
impl<C: alg::RingCommutative> alg::AbstractGroup<alg::Additive> for Cyclotomic<C> {}
impl<C: alg::RingCommutative> alg::AbstractGroupAbelian<alg::Additive> for Cyclotomic<C> {}
// impl<C: alg::RingCommutative> alg::AbstractMonoid<alg::Multiplicative> for Cyclotomic<C> {}
// impl<C: alg::RingCommutative> alg::AbstractRing for Cyclotomic<C> {}
// impl<C: alg::RingCommutative> alg::AbstractRingCommutative for Cyclotomic<C> {}

// impl<C: alg::RingCommutative> alg::RingCommutative for Cyclotomic<C> {}

use alga::general as alg;
use cached::proc_macro::cached;

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

type NN = u32;

#[derive(Clone)]
struct PrimeFactorTerm {
    prime: NN,
    degree: NN,
}

#[cached]
fn factor_impl(prime: NN, degree: NN) -> NN {
    prime.pow(degree as u32)
}

#[cached]
fn factor_phi_impl(prime: NN, degree: NN) -> NN {
    factor_impl(prime, degree) / prime * (prime - 1)
}

impl PrimeFactorTerm {
    fn factor(&self) -> NN {
        factor_impl(self.prime, self.degree)
    }

    fn factor_phi(&self) -> NN {
        factor_phi_impl(self.prime, self.degree)
    }
}

#[cached]
fn factorize(num: NN) -> Vec<PrimeFactorTerm> {
    (2..=num)
        .filter(|&n| (2..n).all(|p| n % p != 0))
        .scan(num, |num, prime| {
            (*num != 1).then(|| PrimeFactorTerm {
                prime,
                degree: std::iter::repeat(prime)
                    .take_while(|prime| (*num % prime == 0).then(|| *num /= prime).is_some())
                    .count() as NN,
            })
        })
        .collect()
    // let mut factorization = Vec::<(usize, usize)>::new();
    // let mut prime = 2usize;
    // while num > 1 {
    //     prime = (prime..)
    //         .find(|&n| factorization.iter().all(|&(p, _)| n % p != 0))
    //         .unwrap();
    //     let mut degree = 0usize;
    //     while num % prime == 0 {
    //         degree += 1;
    //         num /= prime;
    //     }
    //     factorization.push((prime, degree));
    // }
    // factorization
}

#[cached]
fn field_inv(element: NN, order: NN) -> NN {
    // Calculates $denominator^(-1)$ over the finite field of order $order$.
    if element == 1 {
        1
    } else {
        (1 + element * order - field_inv(order % element, element) * order) / element
    }
}

fn field_div(numerator: NN, denominator: NN, order: NN) -> NN {
    // Calculates $numerator * denominator^(-1) % order$ over the finite field of order $order$.
    numerator * field_inv(denominator % order, order) % order
}

// A ring adjoint with roots of unity.
// Basis of coefficients: see Thomas, "Integral Bases for Subfields of Cyclotomic Fields".
// Slightly changed the definition of $J_(k, p)$ to: 1..p if k = 0, 0..p otherwise.
// In the ndarray, each axis represents a factor $p^v$ in $n = prod_i p_i^(v_i)$.
// Elements along an axis represent $zeta_(p^v)^i$ with exponents $i$ running through $p^(v-1)$..$p^v$.
// The size of this axis is thus $p^v - p^(v-1) = phi(p^v)$, and the whole size of ndarray is $phi(n)$.
#[derive(Clone, Debug)]
pub struct Cyclotomic<C> {
    order: NN,
    coeffs: Vec<C>,
}

impl<C> Cyclotomic<C> {
    // pub fn order(&self) -> usize {
    //     self.0.len()
    // }

    pub fn root_of_unity(order: NN, exponent: NN) -> Self
    where
        C: alg::RingCommutative,
    {
        let order_factorization = factorize(order);
        let mut strides = order_factorization
            .iter()
            .rev()
            .scan(1, |stride, &prime_factor_term| {
                *stride *= prime_factor_term.factor_phi();
                Some(*stride)
            })
            .collect::<Vec<_>>();
        strides.rotate_right(1);
        let size = strides
            .get_mut(0)
            .map(|size_ref| std::mem::replace(size_ref, 1))
            .unwrap_or(1);
        strides.reverse();
        let coords = order_factorization
            .iter()
            .map(|&prime_factor_term| {
                let factor = prime_factor_term.factor();
                field_div(exponent, order / factor, factor)
            })
            .collect::<Vec<_>>();
        Self {
            order,
            poly: Polynomial::monomial(exponent),
        }
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

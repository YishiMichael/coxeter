use std::ops::{Range, RangeBounds};

use alga::general as alg;
use cached::proc_macro::cached;
use itertools::Itertools;

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
struct FiniteField {
    prime: NN,
    degree: NN,
    prime_minus_one: NN,
    order_div_prime: NN,
    order: NN,
    order_phi: NN,
}

impl FiniteField {
    fn new(prime: NN, degree: NN) -> Self {
        let prime_minus_one = prime - 1;
        let order_div_prime = prime.pow(degree - 1);
        Self {
            prime,
            degree,
            prime_minus_one,
            order_div_prime,
            order: prime * order_div_prime,
            order_phi: prime_minus_one * order_div_prime,
        }
    }

    fn inv(&self, element: NN) -> NN {
        fn mod_inv(n: NN, m: NN) -> NN {
            if n == 1 {
                1
            } else {
                (1 + n * m - mod_inv(m % n, n) * m) / n
            }
        }
        mod_inv(element, self.order)
    }

    // fn iter_index(&self) -> Range<NN> {
    //     0..self.order_phi
    // }

    // fn index_to_exponent(&self, index: NN) -> NN {
    //     index + self.order_div_prime
    //     //index / self.prime_minus_one * self.prime + index % self.prime_minus_one + 1
    // }

    // fn exponent_to_index_set(&self, exponent: NN) -> (bool, impl Iterator<Item = NN> + '_) {
    //     let flip_sign = exponent < self.order_div_prime;
    //     (
    //         flip_sign,
    //         if flip_sign {
    //             1..=self.prime_minus_one
    //         } else {
    //             0..=0
    //         }
    //         .map(move |k| exponent + k * self.order_div_prime)
    //         .map(|exponent| {
    //             exponent - self.order_div_prime
    //             //exponent / self.prime * self.prime_minus_one + exponent % self.prime - 1
    //         }),
    //     )
    // }
}

#[derive(Clone)]
struct CyclotomicField {
    order: NN,
    finite_fields: Vec<FiniteField>,
    strides: Vec<NN>,
    size: NN,
}

impl CyclotomicField {
    fn new(order: NN) -> Self {
        let finite_fields = (2..=order)
            .filter(|&n| (2..n).all(|p| n % p != 0))
            .scan(order, |num, prime| {
                (*num != 1).then(|| {
                    FiniteField::new(
                        prime,
                        std::iter::repeat(prime)
                            .take_while(|prime| {
                                (*num % prime == 0).then(|| *num /= prime).is_some()
                            })
                            .count() as NN,
                    )
                })
            })
            .collect::<Vec<_>>();
        let mut strides = finite_fields
            .iter()
            .rev()
            .scan(1, |stride, finite_field| {
                *stride *= finite_field.order_phi;
                Some(*stride)
            })
            .collect::<Vec<_>>();
        strides.rotate_right(1);
        let size = strides
            .get_mut(0)
            .map(|size_ref| std::mem::replace(size_ref, 1))
            .unwrap_or(1);
        strides.reverse();
        Self {
            order,
            finite_fields,
            strides,
            size,
        }
    }

    fn decompose_exponent(&self, exponent: NN) -> Vec<NN> {
        self.finite_fields
            .iter()
            .map(|finite_field| {
                exponent * finite_field.inv(self.order / finite_field.order % finite_field.order)
                    % finite_field.order
            })
            .collect::<Vec<_>>()
    }

    // fn index_to_subindices(&self, index: NN) -> impl Iterator<Item = NN> + '_ {
    //     self.strides.iter().scan(index, |state, &stride| {
    //         Some(*state / stride).inspect(|_| *state &= stride)
    //     })
    // }

    fn read_terms<'a, C: alg::RingCommutative + 'a>(
        &'a self,
        coeffs: Vec<C>,
    ) -> impl Iterator<Item = (Vec<NN>, C)> + 'a {
        self.finite_fields
            .iter()
            .map(|finite_field| 0..finite_field.order_phi)
            .multi_cartesian_product()
            .zip(coeffs)
            .filter(|(_, coeff)| !coeff.is_zero())
            .map(move |(indices, coeff)| {
                (
                    self.finite_fields
                        .iter()
                        .zip(indices.into_iter())
                        .map(|(finite_field, index)| index + finite_field.order_div_prime)
                        .collect(),
                    coeff,
                )
            })
    }

    fn write_terms<C: alg::RingCommutative>(
        &self,
        terms: impl Iterator<Item = (Vec<NN>, C)>,
    ) -> Vec<C> {
        let mut coeffs = vec![C::zero(); self.size as usize];
        terms.for_each(|(exponents, coeff)| {
            let mut neg_sign = false;
            let index_sets = self
                .finite_fields
                .iter()
                .zip(exponents.into_iter())
                .map(|(finite_field, exponent)| {
                    let flip_sign = exponent < finite_field.order_div_prime;
                    neg_sign ^= flip_sign;
                    if flip_sign {
                        1..=finite_field.prime_minus_one
                    } else {
                        0..=0
                    }
                    .map(move |k| exponent + k * finite_field.order_div_prime)
                    .map(|exponent| exponent - finite_field.order_div_prime)
                })
                .multi_cartesian_product()
                .collect::<Vec<_>>();
            let signed_coeff = if neg_sign { -coeff } else { coeff };
            index_sets.into_iter().for_each(|indices| {
                let flattened_index = self
                    .strides
                    .iter()
                    .zip(indices.into_iter())
                    .map(|(&stride, index)| index * stride)
                    .sum::<NN>();
                coeffs[flattened_index as usize] = signed_coeff.clone();
            });
        });
        coeffs
    }

    fn embed_coeffs<C: alg::RingCommutative>(
        &self,
        subfield: &'static CyclotomicField,
        coeffs: Vec<C>,
    ) -> Vec<C> {
        self.write_terms(subfield.read_terms(coeffs).map(|(exponents, coeff)| {
            (
                self.finite_fields
                    .iter()
                    .zip_longest(subfield.finite_fields.iter().zip(exponents))
                    .map(|item| match item {
                        itertools::EitherOrBoth::Both(
                            embedded_finite_field,
                            (finite_field, exponent),
                        ) => exponent * (embedded_finite_field.order / finite_field.order),
                        itertools::EitherOrBoth::Left(embedded_finite_field) => {
                            embedded_finite_field.order
                        }
                        _ => unreachable!(),
                    })
                    .collect(),
                coeff,
            )
        }))
        // let mut embedded_coeffs = vec![C::zero(); self.size as usize];
        // cyclotomic_field
        //     .iter_indices()
        //     .zip(coeffs.iter())
        //     .map(|(indices, coeff)| {
        //         self.finite_fields
        //             .iter()
        //             .zip(
        //                 cyclotomic_field
        //                     .finite_fields
        //                     .iter()
        //                     .zip(indices.iter())
        //                     .map(|(field, index)| index + field.order_div_prime),
        //             )
        //             .map(|(embedded_field, exponent)| finite_field.index_to_exponent(index))
        //     })
    }
}

#[cached]
fn cyclotomic_field(order: NN) -> CyclotomicField {
    CyclotomicField::new(order)
}

// fn field_div(numerator: NN, denominator: NN, order: NN) -> NN {
//     // Calculates $numerator * denominator^(-1) % order$ over the finite field of order $order$.
//     numerator * field_inv(denominator % order, order) % order
// }

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

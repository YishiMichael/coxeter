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
    //prime_minus_one: NN,
    order_div_prime: NN,
    order: NN,
    order_phi: NN,
}

impl FiniteField {
    fn new(prime: NN, degree: NN) -> Self {
        //let prime_minus_one = prime - 1;
        let order_div_prime = prime.pow(degree - 1);
        Self {
            prime,
            degree,
            //prime_minus_one,
            order_div_prime,
            order: prime * order_div_prime,
            order_phi: (prime - 1) * order_div_prime,
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
    //strides: Vec<NN>,
    //size: NN,
}

impl CyclotomicField {
    fn new(order: NN) -> Self {
        // let factorization = ;
        // let finite_fields = ;
        // let mut strides = finite_fields
        //     .iter()
        //     .rev()
        //     .scan(1, |stride, finite_field| {
        //         *stride *= finite_field.order_phi;
        //         Some(*stride)
        //     })
        //     .collect::<Vec<_>>();
        // strides.rotate_right(1);
        // let size = strides
        //     .get_mut(0)
        //     .map(|size_ref| std::mem::replace(size_ref, 1))
        //     .unwrap_or(1);
        // strides.reverse();
        Self {
            order,
            finite_fields: std::iter::repeat(())
                .scan((order, 2), |(num, prime), _| {
                    (*num != 1).then(|| (*prime..=*num).find(|factor| *num % *factor == 0).unwrap())
                })
                .dedup_with_count()
                .map(|(degree, prime)| FiniteField::new(prime, degree as NN))
                .collect(),
            //strides,
            //size,
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

    // fn read_terms<'a, C: alg::RingCommutative + 'a>(
    //     &'a self,
    //     coeffs: Vec<C>,
    // ) -> impl Iterator<Item = (Vec<NN>, C)> + 'a {
    //     self.finite_fields
    //         .iter()
    //         .map(|finite_field| 0..finite_field.order_phi)
    //         .multi_cartesian_product()
    //         .zip(coeffs)
    //         .filter(|(_, coeff)| !coeff.is_zero())
    //         .map(move |(indices, coeff)| {
    //             (
    //                 self.finite_fields
    //                     .iter()
    //                     .zip(indices.into_iter())
    //                     .map(|(finite_field, index)| index + finite_field.order_div_prime)
    //                     .collect(),
    //                 coeff,
    //             )
    //         })
    // }

    fn expand_term<C: alg::RingCommutative>(
        &self,
        exponents: Vec<NN>,
        coeff: C,
    ) -> Vec<(Vec<NN>, C)> {
        let mut neg_sign = false;
        let exponents_set = self
            .finite_fields
            .iter()
            .zip(exponents)
            .map(|(finite_field, exponent)| {
                let flip_sign = exponent < finite_field.order_div_prime;
                neg_sign ^= flip_sign;
                if flip_sign {
                    1..finite_field.prime
                } else {
                    0..1
                }
                .map(move |k| exponent + k * finite_field.order_div_prime)
                //.map(|exponent| exponent - finite_field.order_div_prime)
            })
            .multi_cartesian_product()
            .collect::<Vec<_>>();
        let signed_coeff = if neg_sign { -coeff } else { coeff };
        exponents_set
            .into_iter()
            .map(move |exponents| (exponents, signed_coeff.clone()))
            .collect()
    }

    fn embed_terms<C: alg::RingCommutative>(
        &self,
        subfield: &CyclotomicField,
        terms: Vec<(Vec<NN>, C)>,
    ) -> Vec<(Vec<NN>, C)> {
        terms
            .into_iter()
            .map(|(low_exponents, coeff)| {
                let high_exponents = self
                    .finite_fields
                    .iter()
                    .merge_join_by(
                        subfield.finite_fields.iter().zip(low_exponents),
                        |&high_finite_field, &(low_finite_field, _)| {
                            high_finite_field.prime.cmp(&low_finite_field.prime)
                        },
                    )
                    .map(|merge_item| match merge_item {
                        itertools::EitherOrBoth::Both(
                            high_finite_field,
                            (low_finite_field, low_exponent),
                        ) => low_exponent * (high_finite_field.order / low_finite_field.order),
                        itertools::EitherOrBoth::Left(high_finite_field) => high_finite_field.order,
                        _ => unreachable!(),
                    })
                    .collect();
                self.expand_term(high_exponents, coeff)
            })
            .flatten()
            .collect()
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

impl PartialEq for CyclotomicField {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
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

// A ring adjoint with all roots of unity.
// This can be viewed as a graded ring, indexing by factorization exponents.
// Basis of coefficients: see Thomas, "Integral Bases for Subfields of Cyclotomic Fields".
// Slightly changed the definition of $J_(k, p)$ to: 1..p if k = 0, 0..p otherwise.
// In the representation space, each axis represents a factor $p^v$ in $n = prod_i p_i^(v_i)$.
// Elements along an axis represent $zeta_(p^v)^i$ with exponents $i$ running through $p^(v-1)$..$p^v$.
// The size of this axis is thus $p^v - p^(v-1) = phi(p^v)$, and the whole size of representation space is $phi(n)$.
// Invariants: `terms` are arranged in lexigraphical order of exponents; all coefficients are nonzero.
#[derive(Clone, PartialEq)]
pub struct Cyclotomic<C> {
    //order: NN,
    field: CyclotomicField,
    terms: Vec<(Vec<NN>, C)>,
}

impl<C> Cyclotomic<C> {
    // pub fn order(&self) -> usize {
    //     self.0.len()
    // }

    pub fn root_of_unity(order: NN, exponent: NN) -> Self
    where
        C: alg::RingCommutative,
    {
        let field = cyclotomic_field(order);
        let mut terms = field.expand_term(field.decompose_exponent(exponent), C::one());
        terms.sort_by_key(|(exponents, _)| exponents.clone());
        Self { field, terms }
        // let mut number = std::iter::repeat_with(|| <C as alg::Identity<alg::Additive>>::identity()).take(order);
        // number.0[exponent] = <C as alg::Identity<alg::Multiplicative>>::identity();
        // number
    }

    fn embed_terms_to_order(&mut self, order: NN)
    where
        C: alg::RingCommutative,
    {
        if self.field.order % order != 0 {
            let field = cyclotomic_field(num::Integer::lcm(&self.field.order, &order));
            let mut terms = field.embed_terms(&self.field, std::mem::take(&mut self.terms));
            terms.sort_by_key(|(exponents, _)| exponents.clone());
            self.field = field;
            self.terms = terms;
        };
    }
}

// fn zipped<C: alg::RingCommutative>(lhs: Vec<C>, rhs: Vec<C>) -> Vec<(C, C)> {
//     let lhs_order = lhs.len();
//     let rhs_order = rhs.len();
//     let order = num::Integer::lcm(&lhs_order, &rhs_order);
//     let lhs_scalar = order / lhs_order;
//     let rhs_scalar = order / rhs_order;
//     (0..order)
//         .map(|index| {
//             (
//                 if index % lhs_scalar == 0 {
//                     lhs[index / lhs_scalar].clone()
//                 } else {
//                     C::zero()
//                 },
//                 if index % rhs_scalar == 0 {
//                     rhs[index / rhs_scalar].clone()
//                 } else {
//                     C::zero()
//                 },
//             )
//         })
//         .collect()
// }

impl<C: alg::RingCommutative> std::ops::Add for Cyclotomic<C> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<C: alg::RingCommutative> std::ops::AddAssign for Cyclotomic<C> {
    fn add_assign(&mut self, mut rhs: Self) {
        self.embed_terms_to_order(rhs.field.order);
        rhs.embed_terms_to_order(self.field.order);
        let lhs_terms = std::mem::take(&mut self.terms);
        let rhs_terms = rhs.terms;
        self.terms = lhs_terms
            .into_iter()
            .merge_join_by(rhs_terms, |(lhs_exponents, _), (rhs_exponents, _)| {
                lhs_exponents.cmp(rhs_exponents)
            })
            .filter_map(|merge_item| match merge_item {
                itertools::EitherOrBoth::Both((exponents, lhs_coeff), (_, rhs_coeff)) => {
                    let coeff = lhs_coeff + rhs_coeff;
                    (!coeff.is_zero()).then_some((exponents, coeff))
                }
                itertools::EitherOrBoth::Left((exponents, coeff)) => Some((exponents, coeff)),
                itertools::EitherOrBoth::Right((exponents, coeff)) => Some((exponents, coeff)),
            })
            .collect();
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
        *self += -rhs;
    }
}

impl<C: alg::RingCommutative> std::ops::Neg for Cyclotomic<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            field: self.field,
            terms: self
                .terms
                .into_iter()
                .map(|(exponent, coeff)| (exponent, -coeff))
                .collect(),
        }
    }
}

impl<C: alg::RingCommutative> num::Zero for Cyclotomic<C> {
    fn zero() -> Self {
        Self {
            field: cyclotomic_field(1),
            terms: Vec::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
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
        // fn cross_inner_product<C: alg::RingCommutative>(zipped: &[(C, C)]) -> C {
        //     zipped
        //         .iter()
        //         .map(|(elem_0, _)| elem_0)
        //         .zip(zipped.iter().map(|(_, elem_1)| elem_1).rev())
        //         .map(|(elem_0, elem_1)| elem_0.clone() * elem_1.clone())
        //         .fold(C::zero(), C::add)
        // }
        // let zipped = Self::zipped(std::mem::replace(&mut self.0, Vec::new()), rhs.0);
        // self.0.extend((1..=zipped.len()).map(|index| {
        //     cross_inner_product(&zipped[..index]) + cross_inner_product(&zipped[index..])
        // }))
    }
}

impl<C: alg::RingCommutative> num::One for Cyclotomic<C> {
    fn one() -> Self {
        Self::root_of_unity(1, 0)
    }
}

////////////////////////
// impl<C: alg::RingCommutative> PartialEq for Cyclotomic<C> {
//     fn eq(&self, _rhs: &Self) -> bool {
//         todo!()
//     }
// }

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

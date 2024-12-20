use alga::general as alg;
use itertools::Itertools;

trait GradeMeta {
    type Grade: std::cmp::Ord;
}

#[derive(Clone, PartialEq)]
pub struct GradedRing<M, G, C> {
    grade_meta: M,
    terms: Vec<(G, C)>,
}

impl<C: alg::RingCommutative> GradedRing<C> {
    pub fn monomial(exponent: NN) -> Self {
        Self {
            terms: vec![(exponent, C::one())],
        }
    }

    // pub fn simplification(mut self) -> Self
    // where
    //     C: num::Zero,
    // {
    //     while self.coeffs.last().is_some_and(|coeff| coeff.is_zero()) {
    //         self.coeffs.pop().unwrap();
    //     }
    //     self
    // }

    pub fn eval(self, value: C) -> C {
        self.terms
            .into_iter()
            .map(|(exponent, coeff)| {
                std::iter::repeat(value.clone())
                    .take(exponent as usize)
                    .fold(coeff, C::mul)
            })
            .fold(C::zero(), C::add)
    }

    pub fn add_impl(self, rhs: Self) -> Self {}
}

// impl<C: alg::RingCommutative> PartialEq for GradedRing<C> {
//     fn eq(&self, rhs: &Self) -> bool {
//         num::Zero::is_zero(&(self.clone() - rhs.clone()))
//     }
// }

impl<C: alg::RingCommutative> std::ops::Neg for GradedRing<C> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            terms: self
                .terms
                .into_iter()
                .map(|(exponent, coeff)| (exponent, -coeff))
                .collect(),
        }
    }
}

impl<C: alg::RingCommutative> std::ops::AddAssign for GradedRing<C> {
    fn add_assign(&mut self, rhs: Self) {
        let lhs_terms = std::mem::take(&mut self.terms);
        let rhs_terms = rhs.terms;
        self.terms = lhs_terms
            .into_iter()
            .merge_join_by(rhs_terms, |(lhs_exponent, _), (rhs_exponent, _)| {
                lhs_exponent.cmp(rhs_exponent)
            })
            .filter_map(|merge_item| match merge_item {
                itertools::EitherOrBoth::Both((exponent, lhs_coeff), (_, rhs_coeff)) => {
                    let coeff = lhs_coeff + rhs_coeff;
                    (!coeff.is_zero()).then_some((exponent, coeff))
                }
                itertools::EitherOrBoth::Left((exponent, coeff)) => Some((exponent, coeff)),
                itertools::EitherOrBoth::Right((exponent, coeff)) => Some((exponent, coeff)),
            })
            .collect();
    }
}

impl<C: alg::RingCommutative> std::ops::Add for GradedRing<C> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self {
        self += rhs;
        self
    }
}

impl<C: alg::RingCommutative> std::ops::SubAssign for GradedRing<C> {
    fn sub_assign(&mut self, rhs: Self) {
        *self += -rhs;
    }
}

impl<C: alg::RingCommutative> std::ops::Sub for GradedRing<C> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self {
        self -= rhs;
        self
    }
}

impl<C: alg::RingCommutative> num::Zero for GradedRing<C> {
    fn zero() -> Self {
        Self { terms: Vec::new() }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
}

impl<C: alg::RingCommutative> std::ops::MulAssign for GradedRing<C> {
    fn mul_assign(&mut self, rhs: Self) {
        let lhs_terms = std::mem::take(&mut self.terms);
        let rhs_terms = rhs.terms;
        self.terms = lhs_terms
            .into_iter()
            .cartesian_product(rhs_terms)
            .map(|((lhs_exponent, lhs_coeff), (rhs_exponent, rhs_coeff))| ())
            .sorted_by(|(exponent, _)| exponent)
            .collect()
        // (0..self.coeffs.len() + rhs.coeffs.len() - 1).map(|exponent| {
        //     coeffs
        //         .iter()
        //         .enumerate()
        //         .rev()
        //         .filter_map(|(lhs_exponent, lhs_coeff)| {
        //             exponent
        //                 .checked_sub(lhs_exponent)
        //                 .and_then(|rhs_exponent| rhs.coeffs.get(rhs_exponent))
        //                 .map(|rhs_coeff| lhs_coeff.clone() * rhs_coeff.clone())
        //         })
        //         .fold(C::zero(), std::ops::Add::add)
        // });
    }
}

impl<C: alg::RingCommutative> std::ops::Mul for GradedRing<C> {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self {
        self *= rhs;
        self
    }
}

impl<C: alg::RingCommutative> num::One for GradedRing<C> {
    fn one() -> Self {
        Self::monomial(0)
    }
}

impl<C: alg::RingCommutative> alg::AbstractMagma<alg::Additive> for GradedRing<C> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() + rhs.clone()
    }
}
impl<C: alg::RingCommutative> alg::Identity<alg::Additive> for GradedRing<C> {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<C: alg::RingCommutative> alg::TwoSidedInverse<alg::Additive> for GradedRing<C> {
    fn two_sided_inverse(&self) -> Self {
        -self.clone()
    }
}
impl<C: alg::RingCommutative> alg::AbstractMagma<alg::Multiplicative> for GradedRing<C> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() * rhs.clone()
    }
}
impl<C: alg::RingCommutative> alg::Identity<alg::Multiplicative> for GradedRing<C> {
    fn identity() -> Self {
        num::One::one()
    }
}
impl<C: alg::RingCommutative> alg::AbstractSemigroup<alg::Additive> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractMonoid<alg::Additive> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractQuasigroup<alg::Additive> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractLoop<alg::Additive> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractGroup<alg::Additive> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractGroupAbelian<alg::Additive> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractSemigroup<alg::Multiplicative> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractMonoid<alg::Multiplicative> for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractRing for GradedRing<C> {}
impl<C: alg::RingCommutative> alg::AbstractRingCommutative for GradedRing<C> {}

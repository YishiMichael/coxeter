use alga::general as alg;
use itertools::Itertools;

pub trait GradedAlgebraMeta: Clone {
    type Ring: alg::RingCommutative;
    type Grade: alg::AdditiveMonoid + std::cmp::Ord;

    fn annihilate(element: GradedAlgebra<Self>) -> GradedAlgebra<Self> {
        element
    }
}

#[derive(Clone, Debug)]
pub struct GradedAlgebra<M: GradedAlgebraMeta>(Vec<(M::Ring, M::Grade)>);

impl<M: GradedAlgebraMeta> GradedAlgebra<M> {
    pub fn from_terms(terms: impl Iterator<Item = (M::Ring, M::Grade)>) -> Self {
        Self(terms.collect())
    }

    pub fn into_terms(self) -> impl Iterator<Item = (M::Ring, M::Grade)> {
        self.0.into_iter()
    }

    pub fn embed(coeff: M::Ring, grade: M::Grade) -> Self {
        Self::from_terms(std::iter::once((coeff, grade)))
    }

    pub fn embed_scalar(coeff: M::Ring) -> Self {
        Self::embed(coeff, num::Zero::zero())
    }

    pub fn mul_scalar(self, coeff: M::Ring) -> Self {
        self * Self::embed_scalar(coeff)
    }
}

impl<M: GradedAlgebraMeta> PartialEq for GradedAlgebra<M> {
    fn eq(&self, rhs: &Self) -> bool {
        num::Zero::is_zero(&(self.clone() - rhs.clone()))
    }
}

impl<M: GradedAlgebraMeta> std::ops::Neg for GradedAlgebra<M> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::from_terms(self.into_terms().map(|(coeff, grade)| (-coeff, grade)))
    }
}

impl<M: GradedAlgebraMeta> std::ops::Add for GradedAlgebra<M> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::from_terms(
            self.into_terms()
                .merge_join_by(rhs.into_terms(), |(_, grade_0), (_, grade_1)| {
                    grade_0.cmp(grade_1)
                })
                .map(|merge_item| {
                    merge_item
                        .reduce(|(lhs_coeff, grade), (rhs_coeff, _)| (lhs_coeff + rhs_coeff, grade))
                })
                .filter(|(coeff, _)| !num::Zero::is_zero(coeff)),
        )
    }
}

impl<M: GradedAlgebraMeta> std::ops::Sub for GradedAlgebra<M> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl<M: GradedAlgebraMeta> std::ops::Mul for GradedAlgebra<M> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::from_terms(
            self.into_terms()
                .cartesian_product(rhs.into_terms().collect::<Vec<_>>())
                .map(|((lhs_coeff, lhs_grade), (rhs_coeff, rhs_grade))| {
                    (lhs_coeff * rhs_coeff, lhs_grade + rhs_grade)
                })
                .sorted_by_key(|(_, grade)| grade.clone())
                .chunk_by(|(_, grade)| grade.clone())
                .into_iter()
                .map(|(grade, chunk)| {
                    (
                        chunk.fold(num::Zero::zero(), |acc_coeff, (coeff, _)| acc_coeff + coeff),
                        grade.clone(),
                    )
                })
                .filter(|(coeff, _)| !num::Zero::is_zero(coeff)),
        )
    }
}

impl<M: GradedAlgebraMeta> std::ops::AddAssign for GradedAlgebra<M> {
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = lhs + rhs;
    }
}

impl<M: GradedAlgebraMeta> std::ops::SubAssign for GradedAlgebra<M> {
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = lhs - rhs;
    }
}

impl<M: GradedAlgebraMeta> std::ops::MulAssign for GradedAlgebra<M> {
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = lhs * rhs;
    }
}

impl<M: GradedAlgebraMeta> num::Zero for GradedAlgebra<M> {
    fn zero() -> Self {
        Self::from_terms(std::iter::empty())
    }

    fn is_zero(&self) -> bool {
        M::annihilate(self.clone()).into_terms().next().is_none()
    }
}

impl<M: GradedAlgebraMeta> num::One for GradedAlgebra<M> {
    fn one() -> Self {
        Self::embed(num::One::one(), num::Zero::zero())
    }
}

impl<M: GradedAlgebraMeta> alg::AbstractMagma<alg::Additive> for GradedAlgebra<M> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() + rhs.clone()
    }
}
impl<M: GradedAlgebraMeta> alg::Identity<alg::Additive> for GradedAlgebra<M> {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<M: GradedAlgebraMeta> alg::TwoSidedInverse<alg::Additive> for GradedAlgebra<M> {
    fn two_sided_inverse(&self) -> Self {
        -self.clone()
    }
}
impl<M: GradedAlgebraMeta> alg::AbstractMagma<alg::Multiplicative> for GradedAlgebra<M> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() * rhs.clone()
    }
}
impl<M: GradedAlgebraMeta> alg::Identity<alg::Multiplicative> for GradedAlgebra<M> {
    fn identity() -> Self {
        num::One::one()
    }
}
impl<M: GradedAlgebraMeta> alg::AbstractSemigroup<alg::Additive> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractMonoid<alg::Additive> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractQuasigroup<alg::Additive> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractLoop<alg::Additive> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractGroup<alg::Additive> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractGroupAbelian<alg::Additive> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractSemigroup<alg::Multiplicative> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractMonoid<alg::Multiplicative> for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractRing for GradedAlgebra<M> {}
impl<M: GradedAlgebraMeta> alg::AbstractRingCommutative for GradedAlgebra<M> {}

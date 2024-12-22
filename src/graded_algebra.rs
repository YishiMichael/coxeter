use alga::general as alg;
use itertools::Itertools;

pub trait GradedAlgebraMeta: Clone {
    type Ring: alg::RingCommutative;
    type Grade: alg::AdditiveMonoid + std::cmp::Ord;

    fn annihilate(element: GradedAlgebra<Self>) -> GradedAlgebra<Self> {
        element
    }
}

#[derive(Clone)]
pub struct GradedAlgebra<M: GradedAlgebraMeta> {
    terms: Vec<(M::Ring, M::Grade)>,
}

impl<M: GradedAlgebraMeta> GradedAlgebra<M> {
    fn from_terms(terms: impl IntoIterator<Item = (M::Ring, M::Grade)>) -> Self {
        Self {
            terms: terms.into_iter().collect(),
        }
    }

    fn into_terms(self) -> Vec<(M::Ring, M::Grade)> {
        self.terms
    }

    pub fn embed(coeff: M::Ring, grade: M::Grade) -> Self {
        Self::from_terms(std::iter::once((coeff, grade)))
    }
}

impl<M: GradedAlgebraMeta> PartialEq for GradedAlgebra<M> {
    fn eq(&self, other: &Self) -> bool {
        num::Zero::is_zero(&(self.clone() - other.clone()))
    }
}

impl<M: GradedAlgebraMeta> std::iter::Sum for GradedAlgebra<M> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::into_terms)
                .kmerge_by(|(_, grade_0), (_, grade_1)| grade_0 < grade_1)
                .chunk_by(|(_, grade)| grade.clone())
                .into_iter()
                .map(|(grade, chunk)| {
                    (
                        chunk
                            .map(|(coeff, _)| coeff)
                            .fold(<M::Ring as num::Zero>::zero(), |coeff_acc, coeff| {
                                coeff_acc + coeff
                            }),
                        grade.clone(),
                    )
                })
                .filter(|(coeff, _)| !num::Zero::is_zero(coeff)),
        )
    }
}

impl<M: GradedAlgebraMeta> std::iter::Product for GradedAlgebra<M> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::into_terms)
                .multi_cartesian_product()
                .map(|combination| {
                    combination.into_iter().fold(
                        (
                            <M::Ring as num::One>::one(),
                            <M::Grade as num::Zero>::zero(),
                        ),
                        |(coeff_acc, grade_acc), (coeff, grade)| {
                            (coeff_acc * coeff, grade_acc + grade)
                        },
                    )
                })
                .filter(|(coeff, _)| !num::Zero::is_zero(coeff))
                .sorted_by_key(|(_, grade)| grade.clone()),
        )
    }
}

impl<M: GradedAlgebraMeta> std::ops::Neg for GradedAlgebra<M> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::from_terms(
            self.into_terms()
                .into_iter()
                .map(|(coeff, grade)| (-coeff, grade)),
        )
    }
}

impl<M: GradedAlgebraMeta> std::ops::Add for GradedAlgebra<M> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        [self, rhs].into_iter().sum()
    }
}

impl<M: GradedAlgebraMeta> std::ops::AddAssign for GradedAlgebra<M> {
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, <Self as num::Zero>::zero());
        *self = lhs + rhs;
    }
}

impl<M: GradedAlgebraMeta> std::ops::Sub for GradedAlgebra<M> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        [self, -rhs].into_iter().sum()
    }
}

impl<M: GradedAlgebraMeta> std::ops::SubAssign for GradedAlgebra<M> {
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, <Self as num::Zero>::zero());
        *self = lhs - rhs;
    }
}

impl<M: GradedAlgebraMeta> std::ops::Mul for GradedAlgebra<M> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        [self, rhs].into_iter().product()
    }
}

impl<M: GradedAlgebraMeta> std::ops::MulAssign for GradedAlgebra<M> {
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, <Self as num::Zero>::zero());
        *self = lhs * rhs;
    }
}

impl<M: GradedAlgebraMeta> num::Zero for GradedAlgebra<M> {
    fn zero() -> Self {
        Self::from_terms(std::iter::empty())
    }

    fn is_zero(&self) -> bool {
        M::annihilate(self.clone()).into_terms().is_empty()
    }
}

impl<M: GradedAlgebraMeta> num::One for GradedAlgebra<M> {
    fn one() -> Self {
        Self::embed(
            <M::Ring as num::One>::one(),
            <M::Grade as num::Zero>::zero(),
        )
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

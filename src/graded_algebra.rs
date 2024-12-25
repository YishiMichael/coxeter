use alga::general as alg;
use itertools::Itertools;

// pub trait GradedAlgebraMeta: Clone {
//     type Ring: alg::RingCommutative;
//     type Grade: alg::AdditiveMonoid + std::cmp::Ord;

//     fn annihilate(element: GradedAlgebra<Self>) -> GradedAlgebra<Self> {
//         element
//     }
// }

#[derive(Clone)]
pub struct GradedAlgebra<R, G>(Vec<(R, G)>);

impl<R, G> GradedAlgebra<R, G> {
    pub fn from_terms(terms: impl Iterator<Item = (R, G)>) -> Self {
        Self(terms.collect())
    }

    pub fn get_terms(&self) -> impl Iterator<Item = &(R, G)> {
        self.0.iter()
    }

    pub fn take_terms_mut(&mut self) -> impl Iterator<Item = &mut (R, G)> {
        self.0.iter_mut()
    }

    pub fn take_terms(self) -> impl Iterator<Item = (R, G)> {
        self.0.into_iter()
    }

    // pub fn embed(coeff: R, grade: G) -> Self {
    //     Self::from_terms(std::iter::once((coeff, grade)))
    // }

    // pub fn embed_scalar(coeff: R) -> Self {
    //     Self::embed(coeff, num::Zero::zero())
    // }

    // pub fn mul_scalar(self, coeff: R) -> Self {
    //     self * Self::embed_scalar(coeff)
    // }
}

impl<R, G> PartialEq for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn eq(&self, rhs: &Self) -> bool {
        num::Zero::is_zero(&(self.clone() - rhs.clone()))
    }
}

impl<R, G> std::iter::Sum for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::take_terms)
                .kmerge_by(|(_, grade_0), (_, grade_1)| grade_0 < grade_1)
                .chunk_by(|(_, grade)| grade)
                .into_iter()
                .map(|(grade, chunk)| (chunk.map(|(coeff, _)| coeff).sum(), grade))
                // self.into_terms()
                // .merge_join_by(rhs.into_terms(), |(_, grade_0), (_, grade_1)| {
                //     grade_0.cmp(grade_1)
                // })
                // .map(|merge_item| {
                //     merge_item
                //         .reduce(|(lhs_coeff, grade), (rhs_coeff, _)| (lhs_coeff + rhs_coeff, grade))
                // })
                .filter(|(coeff, _)| !num::Zero::is_zero(coeff)),
        )
    }
}

impl<R, G> std::ops::Add for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        [self, rhs].into_iter().sum()
    }
}

impl<R, G> std::ops::AddAssign for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = [lhs, rhs].into_iter().sum()
    }
}

impl<R, G> std::ops::Sub for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        [self, rhs.neg()].into_iter().sum()
    }
}

impl<R, G> std::ops::SubAssign for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = [lhs, rhs.neg()].into_iter().sum()
    }
}

impl<R, G> num::Zero for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn zero() -> Self {
        Self::from_terms(std::iter::empty())
    }

    fn is_zero(&self) -> bool {
        self.get_terms().next().is_none()
    }
}

impl<R, G> std::ops::Neg for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::from_terms(self.into_terms().map(|(coeff, grade)| (coeff.neg(), grade)))
    }
}

impl<R, G> std::iter::Product for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::take_terms)
                .multi_cartesian_product()
                .map(|terms| {
                    let (coeffs, grades) = terms.unzip();
                    (coeffs.product(), grades.sum()) // TODO
                })
                .kmerge_by(|(_, grade_0), (_, grade_1)| grade_0 < grade_1)
                .chunk_by(|(_, grade)| grade)
                .into_iter()
                .map(|(grade, chunk)| (chunk.map(|(coeff, _)| coeff).sum(), grade))
                .filter(|(coeff, _)| !num::Zero::is_zero(coeff)),
            // .sorted_by_key(|(_, grade)| grade.clone())
            // .chunk_by(|(_, grade)| grade.clone())
            // .into_iter()
            // .map(|(grade, chunk)| {
            //     (
            //         chunk.fold(num::Zero::zero(), |acc_coeff, (coeff, _)| acc_coeff + coeff),
            //         grade.clone(),
            //     )
            // })
            // .filter(|(coeff, _)| !num::Zero::is_zero(coeff)),
        )
    }
}

impl<R, G> std::ops::Mul for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        [self, rhs].into_iter().product()
    }
}

impl<R, G> std::ops::MulAssign for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = [lhs, rhs].into_iter().product()
    }
}

impl<R, G> num::One for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn one() -> Self {
        Self::from_terms(std::iter::once((num::One::one(), num::Zero::zero())))
    }

    fn is_one(&self) -> bool {
        self.eq(Self::one())
    }
}

impl<R, G> alg::AbstractMagma<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Add::add(self.clone(), rhs.clone())
    }
}
impl<R, G> alg::Identity<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<R, G> alg::TwoSidedInverse<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn two_sided_inverse(&self) -> Self {
        std::ops::Neg::neg(self.clone())
    }
}
impl<R, G> alg::AbstractMagma<alg::Multiplicative> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Mul::mul(self.clone(), rhs.clone())
    }
}
impl<R, G> alg::Identity<alg::Multiplicative> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
    fn identity() -> Self {
        num::One::one()
    }
}
impl<R, G> alg::AbstractSemigroup<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractMonoid<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractQuasigroup<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractLoop<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractGroup<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractGroupAbelian<alg::Additive> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractSemigroup<alg::Multiplicative> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractMonoid<alg::Multiplicative> for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractRing for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}
impl<R, G> alg::AbstractRingCommutative for GradedAlgebra<R, G>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
    G: alg::AdditiveMonoid + std::iter::Sum + std::cmp::Ord,
{
}

use itertools::Itertools;

use super::alg;
// pub trait GradedAlgebraMeta: Clone {
//     type Ring: alg::RingCommutative;
//     type Grade: alg::AdditiveMonoid + Clone + std::cmp::Ord;

//     fn annihilate(element: GradedAlgebra<Self>) -> GradedAlgebra<Self> {
//         element
//     }
// }

#[derive(Clone, PartialEq)]
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
}

impl<R, G> GradedAlgebra<R, G>
where
    R: alg::Ring + Clone,
    G: alg::AdditiveMonoid + Clone + std::cmp::Ord,
{
    pub fn embed(coeff: R, grade: G) -> Self {
        Self::from_terms(std::iter::once((coeff, grade)))
    }

    pub fn embed_scalar(coeff: R) -> Self {
        Self::embed(coeff, alg::AdditiveIdentity::zero())
    }
}

impl<R, G> alg::AdditiveMagma for GradedAlgebra<R, G>
where
    R: alg::Ring + Clone,
    G: alg::AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::take_terms)
                .kmerge_by(|(_, grade_0), (_, grade_1)| grade_0 < grade_1)
                .chunk_by(|(_, grade)| grade.clone())
                .into_iter()
                .map(|(grade, chunk)| {
                    (
                        alg::AdditiveMagma::sum(chunk.map(|(coeff, _)| coeff)),
                        grade.clone(),
                    )
                })
                // self.into_terms()
                // .merge_join_by(rhs.into_terms(), |(_, grade_0), (_, grade_1)| {
                //     grade_0.cmp(grade_1)
                // })
                // .map(|merge_item| {
                //     merge_item
                //         .reduce(|(lhs_coeff, grade), (rhs_coeff, _)| (lhs_coeff + rhs_coeff, grade))
                // })
                .filter(|(coeff, _)| !coeff.is_zero()),
        )
    }
}

impl<R, G> alg::AdditiveIdentity for GradedAlgebra<R, G>
where
    R: alg::Ring + Clone,
    G: alg::AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn zero() -> Self {
        Self::from_terms(std::iter::empty())
    }

    fn is_zero(&self) -> bool {
        self.get_terms().next().is_none()
    }
}

impl<R, G> alg::AdditiveInverse for GradedAlgebra<R, G>
where
    R: alg::Ring + Clone,
    G: alg::AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn neg(self) -> Self {
        Self::from_terms(self.take_terms().map(|(coeff, grade)| (coeff.neg(), grade)))
    }
}

impl<R, G> alg::MultiplicativeMagma for GradedAlgebra<R, G>
where
    R: alg::Ring + Clone,
    G: alg::AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::take_terms)
                .map(|terms| terms.collect::<Vec<_>>())
                .multi_cartesian_product()
                .map(|terms| {
                    let (coeffs, grades) = terms.into_iter().unzip::<_, _, Vec<_>, Vec<_>>();
                    (
                        alg::MultiplicativeMagma::product(coeffs.into_iter()),
                        alg::AdditiveMagma::sum(grades.into_iter()),
                    )
                })
                .sorted_by_key(|(_, grade)| grade.clone())
                .chunk_by(|(_, grade)| grade.clone())
                .into_iter()
                .map(|(grade, chunk)| {
                    (
                        alg::AdditiveMagma::sum(chunk.map(|(coeff, _)| coeff)),
                        grade.clone(),
                    )
                })
                .filter(|(coeff, _)| !coeff.is_zero()),
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

impl<R, G> alg::MultiplicativeIdentity for GradedAlgebra<R, G>
where
    R: alg::Ring + Clone,
    G: alg::AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn one() -> Self {
        Self::from_terms(std::iter::once((
            alg::MultiplicativeIdentity::one(),
            alg::AdditiveIdentity::zero(),
        )))
    }

    fn is_one(&self) -> bool {
        self.get_terms()
            .exactly_one()
            .is_ok_and(|(coeff, grade)| coeff.is_one() && grade.is_zero())
    }
}

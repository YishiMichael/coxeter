use itertools::Itertools;

use super::alg::*;

#[derive(Clone, Debug, PartialEq)]
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

// impl<R, G> GradedAlgebra<R, G>
// where
//     R: Ring,
//     G: AdditiveMonoid + std::cmp::Ord,
// {
//     pub fn embed_root_of_unity(order: U, exponent: U) -> Self {
//         Self::from_terms(std::iter::once((R::one(), grade)))
//     }
// }

impl<R, G> AdditiveMagma for GradedAlgebra<R, G>
where
    R: Ring + Clone,
    G: AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::take_terms)
                .kmerge_by(|(_, grade_0), (_, grade_1)| grade_0 < grade_1)
                .chunk_by(|(_, grade)| grade.clone())
                .into_iter()
                .map(|(grade, chunk)| (R::sum(chunk.map(|(coeff, _)| coeff)), grade.clone()))
                .filter(|(coeff, _)| !coeff.is_zero()),
        )
    }
}

impl<R, G> AdditiveIdentity for GradedAlgebra<R, G>
where
    R: Ring + Clone,
    G: AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn zero() -> Self {
        Self::from_terms(std::iter::empty())
    }

    fn is_zero(&self) -> bool {
        self.get_terms().next().is_none()
    }
}

impl<R, G> AdditiveInverse for GradedAlgebra<R, G>
where
    R: Ring + Clone,
    G: AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn neg(self) -> Self {
        Self::from_terms(self.take_terms().map(|(coeff, grade)| (coeff.neg(), grade)))
    }
}

impl<R, G> MultiplicativeMagma for GradedAlgebra<R, G>
where
    R: Ring + Clone,
    G: AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn prod<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_terms(
            iter.map(Self::take_terms)
                .map(|terms| terms.collect::<Vec<_>>())
                .multi_cartesian_product()
                .map(|terms| {
                    let (coeffs, grades) = terms.into_iter().unzip::<_, _, Vec<_>, Vec<_>>();
                    (R::prod(coeffs.into_iter()), G::sum(grades.into_iter()))
                })
                .sorted_by_key(|(_, grade)| grade.clone())
                .chunk_by(|(_, grade)| grade.clone())
                .into_iter()
                .map(|(grade, chunk)| (R::sum(chunk.map(|(coeff, _)| coeff)), grade.clone()))
                .filter(|(coeff, _)| !coeff.is_zero()),
        )
    }
}

impl<R, G> MultiplicativeIdentity for GradedAlgebra<R, G>
where
    R: Ring + Clone,
    G: AdditiveMonoid + Clone + std::cmp::Ord,
{
    fn one() -> Self {
        Self::from_terms(std::iter::once((R::one(), G::zero())))
    }

    fn is_one(&self) -> bool {
        self.get_terms()
            .exactly_one()
            .is_ok_and(|(coeff, grade)| coeff.is_one() && grade.is_zero())
    }
}

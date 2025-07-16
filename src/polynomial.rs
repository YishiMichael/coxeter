// use super::alg;
// use super::big_uint::BigUint;
// use super::graded_algebra::GradedAlgebra;

// #[derive(Clone, Debug)]
// pub struct Polynomial<R>(GradedAlgebra<R, BigUint>);

// impl<R> Polynomial<R>
// where
//     R: alg::Ring + Clone,
// {
//     pub fn monomial(coeff: R, exponent: BigUint) -> Self {
//         Self(GradedAlgebra::from_terms(std::iter::once((
//             coeff, exponent,
//         ))))
//     }

//     pub fn eval(&self, value: R) -> R {
//         alg::AdditiveMagma::sum(self.0.get_terms().map(|(coeff, exponent)| {
//             alg::MultiplicativeMagma::product(
//                 std::iter::once(coeff.clone())
//                     .chain(exponent.clone().range().map(|_| value.clone())),
//             )
//         }))
//     }
// }

// impl<R> alg::AdditiveMagma for Polynomial<R>
// where
//     R: alg::Ring + Clone,
// {
//     #[inline]
//     fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(alg::AdditiveMagma::sum(iter.map(|x| x.0)))
//     }
// }

// impl<R> alg::AdditiveIdentity for Polynomial<R>
// where
//     R: alg::Ring + Clone,
// {
//     #[inline]
//     fn zero() -> Self {
//         Self(alg::AdditiveIdentity::zero())
//     }

//     #[inline]
//     fn is_zero(&self) -> bool {
//         self.0.is_zero()
//     }
// }

// impl<R> alg::AdditiveInverse for Polynomial<R>
// where
//     R: alg::Ring + Clone,
// {
//     #[inline]
//     fn neg(self) -> Self {
//         Self(self.0.neg())
//     }
// }

// impl<R> alg::MultiplicativeMagma for Polynomial<R>
// where
//     R: alg::Ring + Clone,
// {
//     #[inline]
//     fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(alg::MultiplicativeMagma::product(iter.map(|x| x.0)))
//     }
// }

// impl<R> alg::MultiplicativeIdentity for Polynomial<R>
// where
//     R: alg::Ring + Clone,
// {
//     #[inline]
//     fn one() -> Self {
//         Self(alg::MultiplicativeIdentity::one())
//     }

//     #[inline]
//     fn is_one(&self) -> bool {
//         self.0.is_one()
//     }
// }

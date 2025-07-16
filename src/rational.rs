// use super::alg;
// // use super::big_uint::BigUint;

// #[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
// pub struct Rational(num::Rational32);

// impl<T> From<T> for Rational
// where
//     T: Into<num::Rational32>,
// {
//     fn from(value: T) -> Self {
//         Self(value.into())
//     }
// }

// impl alg::AdditiveMagma for Rational {
//     #[inline]
//     fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(iter.map(|x| x.0).sum())
//     }
// }

// impl alg::AdditiveIdentity for Rational {
//     #[inline]
//     fn zero() -> Self {
//         Self(num::Zero::zero())
//     }

//     #[inline]
//     fn is_zero(&self) -> bool {
//         num::Zero::is_zero(&self.0)
//     }
// }

// impl alg::AdditiveInverse for Rational {
//     #[inline]
//     fn neg(self) -> Self {
//         Self(std::ops::Neg::neg(self.0))
//     }
// }

// impl alg::MultiplicativeMagma for Rational {
//     #[inline]
//     fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(iter.map(|x| x.0).product())
//     }
// }

// impl alg::MultiplicativeIdentity for Rational {
//     #[inline]
//     fn one() -> Self {
//         Self(num::One::one())
//     }

//     #[inline]
//     fn is_one(&self) -> bool {
//         num::One::is_one(&self.0)
//     }
// }

// impl alg::MultiplicativeInverse for Rational {
//     #[inline]
//     fn inv(self) -> Self {
//         Self(num::traits::Inv::inv(self.0))
//     }
// }

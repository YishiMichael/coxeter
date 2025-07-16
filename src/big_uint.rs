// use super::alg;

// #[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
// pub struct BigUint(num::BigUint);

// impl BigUint {
//     pub fn lcm<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(iter.fold(num::One::one(), |lcm, num| num::Integer::lcm(&lcm, &num.0)))
//     }

//     pub fn range(self) -> BigUintRange {
//         BigUintRange(num::iter::range(num::Zero::zero(), self.0))
//     }

//     pub fn prime_factors(self) -> PrimeFactors {
//         PrimeFactors::from(self)
//     }
// }

// impl From<num::BigUint> for BigUint {
//     fn from(value: num::BigUint) -> Self {
//         Self(value)
//     }
// }

// impl From<usize> for BigUint {
//     fn from(value: usize) -> Self {
//         Self(value.into())
//     }
// }

// impl Into<num::BigUint> for BigUint {
//     fn into(self) -> num::BigUint {
//         self.0
//     }
// }

// pub struct BigUintRange(num::iter::Range<num::BigUint>);

// impl Iterator for BigUintRange {
//     type Item = BigUint;

//     #[inline]
//     fn next(&mut self) -> Option<Self::Item> {
//         self.0.next().map(BigUint)
//     }
// }

// pub struct PrimeFactors {
//     n: num::BigUint,
//     p: num::BigUint,
// }

// impl From<BigUint> for PrimeFactors {
//     fn from(value: BigUint) -> Self {
//         Self {
//             n: value.0,
//             p: num::BigUint::from(2usize),
//         }
//     }
// }

// impl Iterator for PrimeFactors {
//     type Item = BigUint;

//     fn next(&mut self) -> Option<Self::Item> {
//         (!num::One::is_one(&self.n)).then(|| {
//             self.p = num::iter::range_from(self.p.clone())
//                 .find(|p| num::Integer::is_multiple_of(&self.n, p))
//                 .unwrap();
//             self.n /= self.p.clone();
//             BigUint(self.p.clone())
//         })
//     }
// }

// impl alg::AdditiveMagma for BigUint {
//     #[inline]
//     fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(iter.map(|x| x.0).sum())
//     }
// }

// impl alg::AdditiveIdentity for BigUint {
//     #[inline]
//     fn zero() -> Self {
//         Self(num::Zero::zero())
//     }

//     #[inline]
//     fn is_zero(&self) -> bool {
//         num::Zero::is_zero(&self.0)
//     }
// }

// impl alg::MultiplicativeMagma for BigUint {
//     #[inline]
//     fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
//         Self(iter.map(|x| x.0).product())
//     }
// }

// impl alg::MultiplicativeIdentity for BigUint {
//     #[inline]
//     fn one() -> Self {
//         Self(num::One::one())
//     }

//     #[inline]
//     fn is_one(&self) -> bool {
//         num::One::is_one(&self.0)
//     }
// }

use num::Integer;
use num::One;
use num::Zero;

use super::alg;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct BigUint(num::BigUint);

impl BigUint {
    pub fn lcm<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        Self(
            iter.into_iter()
                .fold(num::BigUint::one(), |lcm, num| lcm.lcm(&num.0)),
        )
    }

    pub fn range(self) -> BigUintRange {
        BigUintRange(num::iter::range(num::BigUint::zero(), self.0))
    }

    pub fn prime_factors(self) -> PrimeFactors {
        PrimeFactors {
            n: self.0,
            p: num::BigUint::from(2usize),
        }
    }
}

impl alg::AdditiveMagma for BigUint {
    #[inline]
    fn sum<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        Self(iter.into_iter().map(|x| x.0).sum())
    }
}

impl alg::AdditiveIdentity for BigUint {
    #[inline]
    fn zero() -> Self {
        Self(num::BigUint::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl alg::MultiplicativeMagma for BigUint {
    #[inline]
    fn product<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        Self(iter.into_iter().map(|x| x.0).product())
    }
}

impl alg::MultiplicativeIdentity for BigUint {
    #[inline]
    fn one() -> Self {
        Self(num::BigUint::one())
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

// impl From<num::BigUint> for BigUint {
//     fn from(value: num::BigUint) -> Self {
//         Self(value)
//     }
// }

impl Into<num::BigUint> for BigUint {
    fn into(self) -> num::BigUint {
        self.0
    }
}

pub struct BigUintRange(num::iter::Range<num::BigUint>);

impl Iterator for BigUintRange {
    type Item = BigUint;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(BigUint)
    }
}

pub struct PrimeFactors {
    n: num::BigUint,
    p: num::BigUint,
}

impl Iterator for PrimeFactors {
    type Item = BigUint;

    fn next(&mut self) -> Option<Self::Item> {
        (!self.n.is_one()).then(|| {
            self.p = num::iter::range_from(self.p.clone())
                .find(|p| self.n.is_multiple_of(p))
                .unwrap();
            self.n /= self.p.clone();
            BigUint(self.n.clone())
        })
    }
}

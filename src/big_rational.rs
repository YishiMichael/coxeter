use num::traits::Inv;
use num::One;
use num::Zero;
use std::ops::Neg;

use super::alg;
use super::big_uint::BigUint;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct BigRational(num::BigRational);

impl alg::AdditiveMagma for BigRational {
    #[inline]
    fn sum<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        Self(iter.into_iter().map(|x| x.0).sum())
    }
}

impl alg::AdditiveIdentity for BigRational {
    #[inline]
    fn zero() -> Self {
        Self(num::BigRational::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl alg::AdditiveInverse for BigRational {
    #[inline]
    fn neg(self) -> Self {
        Self(self.0.neg())
    }
}

impl alg::MultiplicativeMagma for BigRational {
    #[inline]
    fn product<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        Self(iter.into_iter().map(|x| x.0).product())
    }
}

impl alg::MultiplicativeIdentity for BigRational {
    #[inline]
    fn one() -> Self {
        Self(num::BigRational::one())
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

impl alg::MultiplicativeInverse for BigRational {
    #[inline]
    fn inv(self) -> Self {
        Self(self.0.inv())
    }
}

impl From<(BigUint, BigUint)> for BigRational {
    fn from((numer, denom): (BigUint, BigUint)) -> Self {
        Self(num::BigRational::new(
            <BigUint as Into<num::BigUint>>::into(numer).into(),
            <BigUint as Into<num::BigUint>>::into(denom).into(),
        ))
    }
}

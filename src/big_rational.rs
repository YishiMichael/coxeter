use super::alg;
use super::big_uint::BigUint;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct BigRational(num::BigRational);

impl From<(BigUint, BigUint)> for BigRational {
    fn from((numer, denom): (BigUint, BigUint)) -> Self {
        Self(num::BigRational::new(
            <BigUint as Into<num::BigUint>>::into(numer).into(),
            <BigUint as Into<num::BigUint>>::into(denom).into(),
        ))
    }
}

impl alg::AdditiveMagma for BigRational {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}

impl alg::AdditiveIdentity for BigRational {
    #[inline]
    fn zero() -> Self {
        Self(num::Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        num::Zero::is_zero(&self.0)
    }
}

impl alg::AdditiveInverse for BigRational {
    #[inline]
    fn neg(self) -> Self {
        Self(std::ops::Neg::neg(self.0))
    }
}

impl alg::MultiplicativeMagma for BigRational {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).product())
    }
}

impl alg::MultiplicativeIdentity for BigRational {
    #[inline]
    fn one() -> Self {
        Self(num::One::one())
    }

    #[inline]
    fn is_one(&self) -> bool {
        num::One::is_one(&self.0)
    }
}

impl alg::MultiplicativeInverse for BigRational {
    #[inline]
    fn inv(self) -> Self {
        Self(num::traits::Inv::inv(self.0))
    }
}

use super::alg;
use super::big_uint::BigUint;

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PruferGroup(num::BigRational);

impl PruferGroup {
    #[inline]
    fn wrap(value: num::BigRational) -> num::BigRational {
        value.clone() - value.floor()
    }
}

impl From<(BigUint, BigUint)> for PruferGroup {
    fn from((numer, denom): (BigUint, BigUint)) -> Self {
        Self(Self::wrap(num::BigRational::new(
            <BigUint as Into<num::BigUint>>::into(numer).into(),
            <BigUint as Into<num::BigUint>>::into(denom).into(),
        )))
    }
}

impl alg::AdditiveMagma for PruferGroup {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(Self::wrap(iter.map(|x| x.0).sum::<num::BigRational>()))
    }
}

impl alg::AdditiveIdentity for PruferGroup {
    #[inline]
    fn zero() -> Self {
        Self(num::Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        num::Zero::is_zero(&self.0)
    }
}

impl alg::AdditiveInverse for PruferGroup {
    #[inline]
    fn neg(self) -> Self {
        Self(Self::wrap(std::ops::Neg::neg(self.0)))
    }
}

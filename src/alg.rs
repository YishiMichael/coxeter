pub trait AdditiveMagma {
    fn sum<I: IntoIterator<Item = Self>>(iter: I) -> Self;
}

pub trait AdditiveIdentity {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
}

pub trait AdditiveInverse {
    fn neg(self) -> Self;
}

pub trait MultiplicativeMagma {
    fn product<I: IntoIterator<Item = Self>>(iter: I) -> Self;
}

pub trait MultiplicativeIdentity {
    fn one() -> Self;
    fn is_one(&self) -> bool;
}

pub trait MultiplicativeInverse {
    fn inv(self) -> Self;
}

pub trait AdditiveMonoid: AdditiveMagma + AdditiveIdentity {}
pub trait AdditiveGroup: AdditiveMonoid + AdditiveInverse {}
pub trait MultiplicativeMonoid: MultiplicativeMagma + MultiplicativeIdentity {}
pub trait MultiplicativeGroup: MultiplicativeMonoid + MultiplicativeInverse {}
pub trait Ring: AdditiveGroup + MultiplicativeMonoid {}
pub trait Field: Ring + MultiplicativeInverse {}

impl<T> AdditiveMonoid for T where T: AdditiveMagma + AdditiveIdentity {}
impl<T> AdditiveGroup for T where T: AdditiveMonoid + AdditiveInverse {}
impl<T> MultiplicativeMonoid for T where T: MultiplicativeMagma + MultiplicativeIdentity {}
impl<T> MultiplicativeGroup for T where T: MultiplicativeMonoid + MultiplicativeInverse {}
impl<T> Ring for T where T: AdditiveGroup + MultiplicativeMonoid {}
impl<T> Field for T where T: Ring + MultiplicativeInverse {}

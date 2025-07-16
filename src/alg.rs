pub trait AdditiveMagma: Sized {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self;

    fn add(self, rhs: Self) -> Self {
        Self::sum([self, rhs].into_iter())
    }

    fn sub(self, rhs: Self) -> Self
    where
        Self: AdditiveInverse,
    {
        self.add(rhs.neg())
    }
}

pub trait AdditiveIdentity {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
}

pub trait AdditiveInverse {
    fn neg(self) -> Self;
}

pub trait MultiplicativeMagma: Sized {
    fn prod<I: Iterator<Item = Self>>(iter: I) -> Self;

    fn mul(self, rhs: Self) -> Self {
        Self::prod([self, rhs].into_iter())
    }

    fn div(self, rhs: Self) -> Self
    where
        Self: MultiplicativeInverse,
    {
        self.mul(rhs.inv())
    }
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

impl<T> AdditiveMagma for T
where
    T: num::Num,
{
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), std::ops::Add::add)
    }
}

impl<T> AdditiveIdentity for T
where
    T: num::Num,
{
    #[inline]
    fn zero() -> Self {
        num::Zero::zero()
    }

    #[inline]
    fn is_zero(&self) -> bool {
        num::Zero::is_zero(self)
    }
}

impl<T> AdditiveInverse for T
where
    T: num::Num,
{
    #[inline]
    fn neg(self) -> Self {
        Self::zero() - self
    }
}

impl<T> MultiplicativeMagma for T
where
    T: num::Num,
{
    #[inline]
    fn prod<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), std::ops::Mul::mul)
    }
}

impl<T> MultiplicativeIdentity for T
where
    T: num::Num,
{
    #[inline]
    fn one() -> Self {
        num::One::one()
    }

    #[inline]
    fn is_one(&self) -> bool {
        num::One::is_one(self)
    }
}

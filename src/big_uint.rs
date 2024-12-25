use alga::general as alg;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct BigUint(num::BigUint);

impl BigUint {
    #[inline]
    pub fn from(value: num::BigUint) -> Self {
        Self(value)
    }

    #[inline]
    pub fn get(&self) -> &num::BigUint {
        &self.0
    }

    #[inline]
    pub fn get_mut(&mut self) -> &mut num::BigUint {
        &mut self.0
    }

    #[inline]
    pub fn take(self) -> num::BigUint {
        self.0
    }
}

impl std::iter::Sum for BigUint {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from(iter.map(|x| x.take()).sum())
    }
}

impl std::ops::Add for BigUint {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::from(self.take().add(rhs.take()))
    }
}

impl std::ops::AddAssign for BigUint {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.get_mut().add_assign(rhs.take());
    }
}

impl num::Zero for BigUint {
    #[inline]
    fn zero() -> Self {
        Self::from(num::BigUint::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.get().is_zero()
    }
}

impl alg::AbstractMagma<alg::Additive> for BigUint {
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Add::add(self.clone(), rhs.clone())
    }
}
impl alg::Identity<alg::Additive> for BigUint {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl alg::AbstractSemigroup<alg::Additive> for BigUint {}
impl alg::AbstractMonoid<alg::Additive> for BigUint {}

use alga::general as alg;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct PruferGroup(num_rational::BigRational);

impl PruferGroup {
    #[inline]
    pub fn from(value: num_rational::BigRational) -> Self {
        Self(value.clone() - value.floor())
    }

    #[inline]
    pub fn get(&self) -> &num_rational::BigRational {
        &self.0
    }

    #[inline]
    pub fn get_mut(&mut self) -> &mut num_rational::BigRational {
        &mut self.0
    }

    #[inline]
    pub fn take(self) -> num_rational::BigRational {
        self.0
    }
}

impl std::iter::Sum for PruferGroup {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from(iter.map(|x| x.take()).sum())
    }
}

impl std::ops::Add for PruferGroup {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::from(self.take().add(rhs.take()))
    }
}

impl std::ops::AddAssign for PruferGroup {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.get_mut().add_assign(rhs.take());
    }
}

impl std::ops::Sub for PruferGroup {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::from(self.take().sub(rhs.take()))
    }
}

impl std::ops::SubAssign for PruferGroup {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.get_mut().sub_assign(rhs.take());
    }
}

impl num::Zero for PruferGroup {
    #[inline]
    fn zero() -> Self {
        Self::from(num_rational::BigRational::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.get().is_zero()
    }
}

impl std::ops::Neg for PruferGroup {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::from(self.take().neg())
    }
}

impl alg::AbstractMagma<alg::Additive> for PruferGroup {
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Add::add(self.clone(), rhs.clone())
    }
}
impl alg::Identity<alg::Additive> for PruferGroup {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl alg::TwoSidedInverse<alg::Additive> for PruferGroup {
    fn two_sided_inverse(&self) -> Self {
        std::ops::Neg::neg(self.clone())
    }
}
impl alg::AbstractSemigroup<alg::Additive> for PruferGroup {}
impl alg::AbstractMonoid<alg::Additive> for PruferGroup {}
impl alg::AbstractQuasigroup<alg::Additive> for PruferGroup {}
impl alg::AbstractLoop<alg::Additive> for PruferGroup {}
impl alg::AbstractGroup<alg::Additive> for PruferGroup {}
impl alg::AbstractGroupAbelian<alg::Additive> for PruferGroup {}

use alga::general as alg;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct BigRational(num_rational::BigRational);

impl BigRational {
    #[inline]
    pub fn from(value: num_rational::BigRational) -> Self {
        Self(value)
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

impl std::iter::Sum for BigRational {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from(iter.map(|x| x.take()).sum())
    }
}

impl std::ops::Add for BigRational {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::from(self.take().add(rhs.take()))
    }
}

impl std::ops::AddAssign for BigRational {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.get_mut().add_assign(rhs.take());
    }
}

impl std::ops::Sub for BigRational {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::from(self.take().sub(rhs.take()))
    }
}

impl std::ops::SubAssign for BigRational {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.get_mut().sub_assign(rhs.take());
    }
}

impl num::Zero for BigRational {
    #[inline]
    fn zero() -> Self {
        Self::from(num_rational::BigRational::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.get().is_zero()
    }
}

impl std::ops::Neg for BigRational {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::from(self.take().neg())
    }
}

impl std::iter::Product for BigRational {
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from(iter.map(|x| x.take()).product())
    }
}

impl std::ops::Mul for BigRational {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::from(self.take().mul(rhs.take()))
    }
}

impl std::ops::MulAssign for BigRational {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.get_mut().mul_assign(rhs.take());
    }
}

impl std::ops::Div for BigRational {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        Self::from(self.take().div(rhs.take()))
    }
}

impl std::ops::DivAssign for BigRational {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.get_mut().div_assign(rhs.take());
    }
}

impl num::One for BigRational {
    #[inline]
    fn one() -> Self {
        Self::from(num_rational::BigRational::one())
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.get().is_one()
    }
}

impl num::traits::Inv for BigRational {
    type Output = Self;

    #[inline]
    fn inv(self) -> Self::Output {
        Self::from(self.take().inv())
    }
}

impl alg::AbstractMagma<alg::Additive> for BigRational {
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Add::add(self.clone(), rhs.clone())
    }
}
impl alg::Identity<alg::Additive> for BigRational {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl alg::TwoSidedInverse<alg::Additive> for BigRational {
    fn two_sided_inverse(&self) -> Self {
        std::ops::Neg::neg(self.clone())
    }
}
impl alg::AbstractMagma<alg::Multiplicative> for BigRational {
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Mul::mul(self.clone(), rhs.clone())
    }
}
impl alg::Identity<alg::Multiplicative> for BigRational {
    fn identity() -> Self {
        num::One::one()
    }
}
impl alg::TwoSidedInverse<alg::Multiplicative> for BigRational {
    fn two_sided_inverse(&self) -> Self {
        num::traits::Inv::inv(self.clone())
    }
}
impl alg::AbstractSemigroup<alg::Additive> for BigRational {}
impl alg::AbstractMonoid<alg::Additive> for BigRational {}
impl alg::AbstractQuasigroup<alg::Additive> for BigRational {}
impl alg::AbstractLoop<alg::Additive> for BigRational {}
impl alg::AbstractGroup<alg::Additive> for BigRational {}
impl alg::AbstractGroupAbelian<alg::Additive> for BigRational {}
impl alg::AbstractSemigroup<alg::Multiplicative> for BigRational {}
impl alg::AbstractMonoid<alg::Multiplicative> for BigRational {}
impl alg::AbstractQuasigroup<alg::Multiplicative> for BigRational {}
impl alg::AbstractLoop<alg::Multiplicative> for BigRational {}
impl alg::AbstractGroup<alg::Multiplicative> for BigRational {}
impl alg::AbstractGroupAbelian<alg::Multiplicative> for BigRational {}
impl alg::AbstractRing for BigRational {}
impl alg::AbstractRingCommutative for BigRational {}
impl alg::AbstractField for BigRational {}

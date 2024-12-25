use alga::general as alg;

#[derive(Clone)]
pub struct BigRational(num_rational::BigRational);

impl BigRational {
    pub fn new(numer: num::BigInt, denom: num::BigInt) -> Self {
        Self(num_rational::BigRational::new(numer, denom))
    }
}

impl PartialEq for BigRational {
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl Eq for BigRational {}

impl PartialOrd for BigRational {
    fn partial_cmp(&self, rhs: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&rhs.0)
    }
}

impl Ord for BigRational {
    fn cmp(&self, rhs: &Self) -> std::cmp::Ordering {
        self.0.cmp(&rhs.0)
    }
}

impl std::iter::Sum for BigRational {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).sum())
    }
}

impl std::ops::Add for BigRational {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0.add(rhs.0))
    }
}

impl std::ops::AddAssign for BigRational {
    fn add_assign(&mut self, rhs: Self) {
        self.0.add_assign(rhs.0);
    }
}

impl std::ops::Sub for BigRational {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0.sub(rhs.0))
    }
}

impl std::ops::SubAssign for BigRational {
    fn sub_assign(&mut self, rhs: Self) {
        self.0.sub_assign(rhs.0);
    }
}

impl num::Zero for BigRational {
    fn zero() -> Self {
        Self(num_rational::BigRational::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    fn set_zero(&mut self) {
        self.0.set_zero();
    }
}

impl std::ops::Neg for BigRational {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl std::iter::Product for BigRational {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self(iter.map(|x| x.0).product())
    }
}

impl std::ops::Mul for BigRational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.mul(rhs.0))
    }
}

impl std::ops::MulAssign for BigRational {
    fn mul_assign(&mut self, rhs: Self) {
        self.0.mul_assign(rhs.0);
    }
}

impl std::ops::Div for BigRational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0.div(rhs.0))
    }
}

impl std::ops::DivAssign for BigRational {
    fn div_assign(&mut self, rhs: Self) {
        self.0.div_assign(rhs.0);
    }
}

impl num::One for BigRational {
    fn one() -> Self {
        Self(num_rational::BigRational::one())
    }

    fn is_one(&self) -> bool {
        self.0.is_one()
    }

    fn set_one(&mut self) {
        self.0.set_one();
    }
}

impl num::traits::Inv for BigRational {
    type Output = Self;

    fn inv(self) -> Self::Output {
        Self(self.0.inv())
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

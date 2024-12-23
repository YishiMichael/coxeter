use alga::general as alg;
use num::rational::Ratio;

#[derive(Clone, Copy, Debug)]
pub struct PruferMonoid<T>(Ratio<T>);

impl<T: Clone + num::Integer> PruferMonoid<T> {
    #[inline]
    pub fn new(ratio: Ratio<T>) -> Self {
        Self(ratio.clone() - ratio.floor())
    }

    #[inline]
    pub const fn numer(&self) -> &T {
        self.0.numer()
    }

    #[inline]
    pub const fn denom(&self) -> &T {
        self.0.denom()
    }
}

impl<T: Clone + num::Integer> Ord for PruferMonoid<T> {
    #[inline]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

impl<T: Clone + num::Integer> PartialOrd for PruferMonoid<T> {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<T: Clone + num::Integer> PartialEq for PruferMonoid<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl<T: Clone + num::Integer> Eq for PruferMonoid<T> {}

impl<T: Clone + num::Integer> std::ops::Add for PruferMonoid<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.0.add(rhs.0))
    }
}

impl<T: Clone + num::Integer> std::ops::Sub for PruferMonoid<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.0.sub(rhs.0))
    }
}

impl<T: Clone + num::Integer + num::traits::NumAssign> std::ops::AddAssign for PruferMonoid<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<T: Clone + num::Integer + num::traits::NumAssign> std::ops::SubAssign for PruferMonoid<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<T: Clone + num::Integer> num::Zero for PruferMonoid<T> {
    #[inline]
    fn zero() -> Self {
        Self::new(Ratio::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[inline]
    fn set_zero(&mut self) {
        self.0.set_zero();
    }
}

impl<T: Clone + num::Integer> alg::AbstractMagma<alg::Additive> for PruferMonoid<T> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() + rhs.clone()
    }
}
impl<T: Clone + num::Integer> alg::Identity<alg::Additive> for PruferMonoid<T> {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<T: Clone + num::Integer> alg::AbstractSemigroup<alg::Additive> for PruferMonoid<T> {}
impl<T: Clone + num::Integer> alg::AbstractMonoid<alg::Additive> for PruferMonoid<T> {}

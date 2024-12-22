use alga::general as alg;

// natural number
pub trait NN: num::Integer + num::Unsigned + Clone + Copy + std::cmp::Ord {}

impl NN for u32 {}

// Always be a reduced fraction. Zero is represented as $0 / 1$.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PruferGroup<N: NN> {
    pub denominator: N,
    pub numerator: N,
}

impl<N: NN> PruferGroup<N> {
    pub fn new(denominator: N, numerator: N) -> Self {
        if numerator.is_zero() {
            Self {
                denominator: N::one(),
                numerator: N::zero(),
            }
        } else {
            let numerator = numerator % denominator;
            let gcd = numerator.gcd(&denominator);
            Self {
                denominator: denominator / gcd,
                numerator: numerator / gcd,
            }
        }
    }
}

impl<N: NN> std::ops::Neg for PruferGroup<N> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(self.denominator, self.denominator - self.numerator)
    }
}

impl<N: NN> std::ops::Add for PruferGroup<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let (gcd, denominator) = self.denominator.gcd_lcm(&rhs.denominator);
        let numerator = (self.numerator * rhs.denominator + rhs.numerator * self.denominator) / gcd;
        Self::new(denominator, numerator)
    }
}

impl<N: NN> std::ops::Sub for PruferGroup<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl<N: NN> std::ops::AddAssign for PruferGroup<N> {
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, <Self as num::Zero>::zero());
        *self = lhs + rhs;
    }
}

impl<N: NN> std::ops::SubAssign for PruferGroup<N> {
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, <Self as num::Zero>::zero());
        *self = lhs - rhs;
    }
}

impl<N: NN> num::Zero for PruferGroup<N> {
    fn zero() -> Self {
        Self::new(N::one(), N::zero())
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

impl<N: NN> alg::AbstractMagma<alg::Additive> for PruferGroup<N> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() + rhs.clone()
    }
}
impl<N: NN> alg::Identity<alg::Additive> for PruferGroup<N> {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<N: NN> alg::TwoSidedInverse<alg::Additive> for PruferGroup<N> {
    fn two_sided_inverse(&self) -> Self {
        -self.clone()
    }
}
impl<N: NN> alg::AbstractSemigroup<alg::Additive> for PruferGroup<N> {}
impl<N: NN> alg::AbstractMonoid<alg::Additive> for PruferGroup<N> {}
impl<N: NN> alg::AbstractQuasigroup<alg::Additive> for PruferGroup<N> {}
impl<N: NN> alg::AbstractLoop<alg::Additive> for PruferGroup<N> {}
impl<N: NN> alg::AbstractGroup<alg::Additive> for PruferGroup<N> {}
impl<N: NN> alg::AbstractGroupAbelian<alg::Additive> for PruferGroup<N> {}

use alga::general as alg;
use itertools::Itertools;

#[derive(Clone, Debug)]
pub struct Polynomial<C> {
    coeffs: Vec<C>,
}

impl<C: alg::RingCommutative> Polynomial<C> {
    pub fn monomial(exponent: usize) -> Self {
        Self {
            coeffs: std::iter::repeat_with(|| C::zero())
                .take(exponent)
                .chain(std::iter::once(C::one()))
                .collect(),
        }
    }

    // pub fn simplification(mut self) -> Self
    // where
    //     C: num::Zero,
    // {
    //     while self.coeffs.last().is_some_and(|coeff| coeff.is_zero()) {
    //         self.coeffs.pop().unwrap();
    //     }
    //     self
    // }

    pub fn eval(self, value: C) -> C {
        self.coeffs
            .into_iter()
            .rev()
            .fold(C::zero(), |acc, coeff| acc * value.clone() + coeff)
    }
}

impl<C: alg::RingCommutative> PartialEq for Polynomial<C> {
    fn eq(&self, rhs: &Self) -> bool {
        num::Zero::is_zero(&(self.clone() - rhs.clone()))
    }
}

impl<C: alg::RingCommutative> std::ops::Neg for Polynomial<C> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            coeffs: self
                .coeffs
                .iter()
                .cloned()
                .map(std::ops::Neg::neg)
                .collect(),
        }
    }
}

impl<C: alg::RingCommutative> std::ops::AddAssign for Polynomial<C> {
    fn add_assign(&mut self, rhs: Self) {
        let coeffs = std::mem::replace(&mut self.coeffs, Vec::new());
        self.coeffs.extend(
            coeffs
                .into_iter()
                .zip_longest(rhs.coeffs.into_iter())
                .map(|pair| pair.reduce(std::ops::Add::add)),
        );
    }
}

impl<C: alg::RingCommutative> std::ops::Add for Polynomial<C> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self {
        self += rhs;
        self
    }
}

impl<C: alg::RingCommutative> std::ops::SubAssign for Polynomial<C> {
    fn sub_assign(&mut self, rhs: Self) {
        *self += -rhs;
    }
}

impl<C: alg::RingCommutative> std::ops::Sub for Polynomial<C> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self {
        self -= rhs;
        self
    }
}

impl<C: alg::RingCommutative> num::Zero for Polynomial<C> {
    fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(num::Zero::is_zero)
    }
}

impl<C: alg::RingCommutative> std::ops::MulAssign for Polynomial<C> {
    fn mul_assign(&mut self, rhs: Self) {
        let coeffs = std::mem::replace(&mut self.coeffs, Vec::new());
        self.coeffs.extend(
            (0..self.coeffs.len() + rhs.coeffs.len() - 1).map(|exponent| {
                coeffs
                    .iter()
                    .enumerate()
                    .rev()
                    .filter_map(|(lhs_exponent, lhs_coeff)| {
                        exponent
                            .checked_sub(lhs_exponent)
                            .and_then(|rhs_exponent| rhs.coeffs.get(rhs_exponent))
                            .map(|rhs_coeff| lhs_coeff.clone() * rhs_coeff.clone())
                    })
                    .fold(C::zero(), std::ops::Add::add)
            }),
        );
    }
}

impl<C: alg::RingCommutative> std::ops::Mul for Polynomial<C> {
    type Output = Self;
    fn mul(mut self, rhs: Self) -> Self {
        self *= rhs;
        self
    }
}

impl<C: alg::RingCommutative> num::One for Polynomial<C> {
    fn one() -> Self {
        Self::monomial(1)
    }
}

impl<C: alg::RingCommutative> alg::AbstractMagma<alg::Additive> for Polynomial<C> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() + rhs.clone()
    }
}
impl<C: alg::RingCommutative> alg::Identity<alg::Additive> for Polynomial<C> {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<C: alg::RingCommutative> alg::TwoSidedInverse<alg::Additive> for Polynomial<C> {
    fn two_sided_inverse(&self) -> Self {
        -self.clone()
    }
}
impl<C: alg::RingCommutative> alg::AbstractMagma<alg::Multiplicative> for Polynomial<C> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() * rhs.clone()
    }
}
impl<C: alg::RingCommutative> alg::Identity<alg::Multiplicative> for Polynomial<C> {
    fn identity() -> Self {
        num::One::one()
    }
}
impl<C: alg::RingCommutative> alg::AbstractSemigroup<alg::Additive> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractMonoid<alg::Additive> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractQuasigroup<alg::Additive> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractLoop<alg::Additive> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractGroup<alg::Additive> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractGroupAbelian<alg::Additive> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractSemigroup<alg::Multiplicative> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractMonoid<alg::Multiplicative> for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractRing for Polynomial<C> {}
impl<C: alg::RingCommutative> alg::AbstractRingCommutative for Polynomial<C> {}

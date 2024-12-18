use alga::general as alg;
use itertools::Itertools;

#[derive(Clone, Debug)]
pub struct Polynomial<C>(Vec<C>);

impl<C> Polynomial<C> {
    pub fn monomial(degree: usize) -> Self
    where
        C: alg::Identity<alg::Additive> + alg::Identity<alg::Multiplicative>,
    {
        Self(
            std::iter::repeat_with(|| <C as alg::Identity<alg::Additive>>::identity())
                .take(degree)
                .chain(std::iter::once(
                    <C as alg::Identity<alg::Multiplicative>>::identity(),
                ))
                .collect(),
        )
    }

    pub fn simplification(mut self) -> Self
    where
        C: num::Zero,
    {
        while self.0.last().is_some_and(|coeff| coeff.is_zero()) {
            self.0.pop().unwrap();
        }
        self
    }

    pub fn eval(self, value: C) -> C
    where
        C: alg::RingCommutative,
    {
        self.0
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
        Self(self.0.iter().cloned().map(std::ops::Neg::neg).collect())
    }
}

impl<C: alg::RingCommutative> std::ops::AddAssign for Polynomial<C> {
    fn add_assign(&mut self, rhs: Self) {
        let coeffs = std::mem::replace(&mut self.0, Vec::new());
        self.0.extend(
            coeffs
                .into_iter()
                .zip_longest(rhs.0.into_iter())
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
        Self(Vec::new())
    }

    fn is_zero(&self) -> bool {
        self.clone().simplification().0.is_empty()
    }
}

impl<C: alg::RingCommutative> std::ops::MulAssign for Polynomial<C> {
    fn mul_assign(&mut self, rhs: Self) {
        let coeffs = std::mem::replace(&mut self.0, Vec::new());
        self.0
            .extend((0..self.0.len() + rhs.0.len() - 1).map(|exp| {
                coeffs
                    .iter()
                    .enumerate()
                    .rev()
                    .filter_map(|(lhs_exp, lhs_coeff)| {
                        exp.checked_sub(lhs_exp)
                            .and_then(|rhs_exp| rhs.0.get(rhs_exp))
                            .map(|rhs_coeff| lhs_coeff.clone() * rhs_coeff.clone())
                    })
                    .fold(C::zero(), std::ops::Add::add)
            }));
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

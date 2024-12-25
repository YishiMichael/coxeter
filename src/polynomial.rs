use alga::general as alg;

use super::big_uint::BigUint;
use super::graded_algebra::GradedAlgebra;

#[derive(Clone, PartialEq)]
pub struct Polynomial<R>(GradedAlgebra<R, BigUint>);

impl<R> Polynomial<R> {
    #[inline]
    pub fn from(value: GradedAlgebra<R, BigUint>) -> Self {
        Self(value)
    }

    #[inline]
    pub fn get(&self) -> &GradedAlgebra<R, BigUint> {
        &self.0
    }

    #[inline]
    pub fn get_mut(&mut self) -> &mut GradedAlgebra<R, BigUint> {
        &mut self.0
    }

    #[inline]
    pub fn take(self) -> GradedAlgebra<R, BigUint> {
        self.0
    }
}

impl<R> Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    pub fn monomial(coeff: R, exponent: BigUint) -> Self {
        Self::from(GradedAlgebra::from_terms(std::iter::once((
            coeff, exponent,
        ))))
    }

    pub fn eval(&self, value: R) -> R {
        self.get()
            .get_terms()
            .map(|(coeff, exponent)| {
                std::iter::repeat(value.clone())
                    .take(exponent.into())
                    .fold(coeff, R::mul)
            })
            .fold(R::zero(), R::add)
    }
}

impl<R> std::iter::Sum for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from(iter.map(|x| x.take()).sum())
    }
}

impl<R> std::ops::Add for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::from(self.take().add(rhs.take()))
    }
}

impl<R> std::ops::AddAssign for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.get_mut().add_assign(rhs.take());
    }
}

impl<R> std::ops::Sub for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::from(self.take().sub(rhs.take()))
    }
}

impl<R> std::ops::SubAssign for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.get_mut().sub_assign(rhs.take());
    }
}

impl<R> num::Zero for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn zero() -> Self {
        Self::from(num_rational::BigRational::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.get().is_zero()
    }
}

impl<R> std::ops::Neg for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::from(self.take().neg())
    }
}

impl<R> std::iter::Product for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from(iter.map(|x| x.take()).product())
    }
}

impl<R> std::ops::Mul for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::from(self.take().mul(rhs.take()))
    }
}

impl<R> std::ops::MulAssign for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.get_mut().mul_assign(rhs.take());
    }
}

impl<R> num::One for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    #[inline]
    fn one() -> Self {
        Self::from(num_rational::BigRational::one())
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.get().is_one()
    }
}

impl<R> alg::AbstractMagma<alg::Additive> for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Add::add(self.clone(), rhs.clone())
    }
}
impl<R> alg::Identity<alg::Additive> for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<R> alg::TwoSidedInverse<alg::Additive> for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    fn two_sided_inverse(&self) -> Self {
        std::ops::Neg::neg(self.clone())
    }
}
impl<R> alg::AbstractMagma<alg::Multiplicative> for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    fn operate(&self, rhs: &Self) -> Self {
        std::ops::Mul::mul(self.clone(), rhs.clone())
    }
}
impl<R> alg::Identity<alg::Multiplicative> for Polynomial<R>
where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product,
{
    fn identity() -> Self {
        num::One::one()
    }
}
impl<R> alg::AbstractSemigroup<alg::Additive> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractMonoid<alg::Additive> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractQuasigroup<alg::Additive> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractLoop<alg::Additive> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractGroup<alg::Additive> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractGroupAbelian<alg::Additive> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractSemigroup<alg::Multiplicative> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractMonoid<alg::Multiplicative> for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractRing for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}
impl<R> alg::AbstractRingCommutative for Polynomial<R> where
    R: alg::RingCommutative + std::iter::Sum + std::iter::Product
{
}

// #[derive(Clone, Debug)]
// pub struct PolynomialMeta<R> {
//     phantom: std::marker::PhantomData<R>,
// }

// impl<R> GradedAlgebraMeta for PolynomialMeta<R>
// where
//     R: alg::RingCommutative,
// {
//     type Ring = R;
//     type Grade = u64;
// }

// pub type Polynomial<R> = GradedAlgebra<PolynomialMeta<R>>;

// impl<R: alg::RingCommutative> Polynomial<R> {
//     pub fn monomial(coeff: R, exponent: u64) -> Self {
//         Self::embed(coeff, exponent)
//     }

//     pub fn eval(self, value: R) -> R {
//         self.into_terms()
//             .map(|(coeff, exponent)| {
//                 std::iter::repeat(value.clone())
//                     .take(exponent as usize)
//                     .fold(coeff, R::mul)
//             })
//             .fold(R::zero(), R::add)
//     }
// }

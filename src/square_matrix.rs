use alga::general as alg;
use itertools::Itertools;
use petgraph::visit::NodeIndexable;

use super::polynomial::Polynomial;

#[derive(Clone)]
pub struct SquareMatrix<T, const N: usize>(
    petgraph::matrix_graph::MatrixGraph<(), T, petgraph::Directed, Option<T>, usize>,
);

impl<T: alg::Ring, const N: usize> SquareMatrix<T, N> {
    pub fn from_fn<F>(mut f: F) -> Self
    where
        F: FnMut(usize, usize) -> T,
    {
        Self(petgraph::matrix_graph::MatrixGraph::from_edges(
            (0..N).cartesian_product(0..N).map(|(i, j)| (i, j, f(i, j))),
        ))
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        self.0
            .edge_weight(self.0.from_index(i), self.0.from_index(j))
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        self.0
            .edge_weight_mut(self.0.from_index(i), self.0.from_index(j))
    }

    pub fn replace(&mut self, i: usize, j: usize, src: T) -> T {
        std::mem::replace(self.get_mut(i, j), src)
    }

    pub fn take(&mut self, i: usize, j: usize) -> T {
        self.replace(i, j, T::zero())
    }

    pub fn transpose(mut self) -> Self {
        Self::from_fn(|i, j| self.take(j, i))
    }

    pub fn determinant(&self) -> T {
        (0..N)
            .permutations(N)
            .map(|permutation| {
                let neg_sign = permutation
                    .iter()
                    .combinations(2)
                    .map(|index_pair| {
                        if let [i, j] = index_pair.as_slice() {
                            i > j
                        } else {
                            unreachable!();
                        }
                    })
                    .fold(false, std::ops::BitXor::bitxor);
                let product = (0..N)
                    .map(|i| self.get(i, permutation[i]))
                    .fold(T::one(), |product, entry| product * entry.clone());
                if neg_sign {
                    -product
                } else {
                    product
                }
            })
            .fold(T::zero(), |product_acc, product| product_acc + product)
    }
}

impl<T: alg::RingCommutative, const N: usize> SquareMatrix<T, N> {
    pub fn characteristic_polynomial(&self) -> Polynomial<T> {
        SquareMatrix::<Polynomial<T>, N>::from_fn(|i, j| {
            (if i == j {
                Polynomial::monomial(T::one(), 1)
            } else {
                num::Zero::zero()
            }) - Polynomial::monomial(self.get(i, j).clone(), 0)
        })
        .determinant()
    }
}

impl<T: alg::Ring, const N: usize> PartialEq for SquareMatrix<T, N> {
    fn eq(&self, rhs: &Self) -> bool {
        (0..N)
            .cartesian_product(0..N)
            .all(|(i, j)| self.get(i, j) == rhs.get(i, j))
    }
}

impl<T: alg::Ring, const N: usize> std::ops::Neg for SquareMatrix<T, N> {
    type Output = Self;

    fn neg(mut self) -> Self {
        Self::from_fn(|i, j| -self.take(i, j))
    }
}

impl<T: alg::Ring, const N: usize> std::ops::Add for SquareMatrix<T, N> {
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self {
        Self::from_fn(|i, j| self.take(i, j) + rhs.take(i, j))
    }
}

impl<T: alg::Ring, const N: usize> std::ops::Sub for SquareMatrix<T, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl<T: alg::Ring, const N: usize> std::ops::Mul for SquareMatrix<T, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::from_fn(|i, j| {
            (0..N)
                .map(|k| self.get(i, k).clone() * rhs.get(k, j).clone())
                .fold(T::zero(), T::add)
        })
    }
}

impl<T: alg::Ring, const N: usize> std::ops::AddAssign for SquareMatrix<T, N> {
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = lhs + rhs;
    }
}

impl<T: alg::Ring, const N: usize> std::ops::SubAssign for SquareMatrix<T, N> {
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = lhs - rhs;
    }
}

impl<T: alg::Ring, const N: usize> std::ops::MulAssign for SquareMatrix<T, N> {
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, num::Zero::zero());
        *self = lhs * rhs;
    }
}

impl<T: alg::Ring, const N: usize> num::Zero for SquareMatrix<T, N> {
    fn zero() -> Self {
        Self::from_fn(|_, _| T::zero())
    }

    fn is_zero(&self) -> bool {
        (0..N)
            .cartesian_product(0..N)
            .all(|(i, j)| self.get(i, j).is_zero())
    }
}

impl<T: alg::Ring, const N: usize> num::One for SquareMatrix<T, N> {
    fn one() -> Self {
        Self::from_fn(|i, j| if i == j { T::one() } else { T::zero() })
    }
}

impl<T: alg::Ring, const N: usize> alg::AbstractMagma<alg::Additive> for SquareMatrix<T, N> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() + rhs.clone()
    }
}
impl<T: alg::Ring, const N: usize> alg::Identity<alg::Additive> for SquareMatrix<T, N> {
    fn identity() -> Self {
        num::Zero::zero()
    }
}
impl<T: alg::Ring, const N: usize> alg::TwoSidedInverse<alg::Additive> for SquareMatrix<T, N> {
    fn two_sided_inverse(&self) -> Self {
        -self.clone()
    }
}
impl<T: alg::Ring, const N: usize> alg::AbstractMagma<alg::Multiplicative> for SquareMatrix<T, N> {
    fn operate(&self, rhs: &Self) -> Self {
        self.clone() * rhs.clone()
    }
}
impl<T: alg::Ring, const N: usize> alg::Identity<alg::Multiplicative> for SquareMatrix<T, N> {
    fn identity() -> Self {
        num::One::one()
    }
}
impl<T: alg::Ring, const N: usize> alg::AbstractSemigroup<alg::Additive> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractMonoid<alg::Additive> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractQuasigroup<alg::Additive> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractLoop<alg::Additive> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractGroup<alg::Additive> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractGroupAbelian<alg::Additive> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractSemigroup<alg::Multiplicative>
    for SquareMatrix<T, N>
{
}
impl<T: alg::Ring, const N: usize> alg::AbstractMonoid<alg::Multiplicative> for SquareMatrix<T, N> {}
impl<T: alg::Ring, const N: usize> alg::AbstractRing for SquareMatrix<T, N> {}

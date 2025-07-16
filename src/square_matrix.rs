use std::fmt::Debug;

use itertools::Itertools;

use super::alg::*;
// use super::polynomial::Polynomial;

#[derive(Clone, Debug)]
pub struct SquareMatrix<T> {
    dimension: usize,
    elements: Vec<Vec<T>>,
}

impl<T> SquareMatrix<T>
where
    T: Ring + Clone,
{
    pub fn from_fn<F>(dimension: usize, mut f: F) -> Self
    where
        F: FnMut(usize, usize) -> T,
    {
        Self {
            dimension,
            elements: (0..dimension)
                .map(|i| (0..dimension).map(|j| f(i, j)).collect())
                .collect(),
        }
    }

    pub fn embed(dimension: usize, value: T) -> Self {
        Self::from_fn(
            dimension,
            |i, j| {
                if i == j {
                    value.clone()
                } else {
                    T::zero()
                }
            },
        )
    }

    pub fn zero(dimension: usize) -> Self {
        Self::from_fn(dimension, |_, _| T::zero())
    }

    pub fn is_zero(&self) -> bool {
        (0..self.dimension())
            .cartesian_product(0..self.dimension())
            .all(|(i, j)| self.get(i, j).is_zero())
    }

    pub fn one(dimension: usize) -> Self {
        Self::from_fn(dimension, |i, j| if i == j { T::one() } else { T::zero() })
    }

    pub fn is_one(&self) -> bool {
        (0..self.dimension())
            .cartesian_product(0..self.dimension())
            .all(|(i, j)| {
                if i == j {
                    self.get(i, j).is_one()
                } else {
                    self.get(i, j).is_zero()
                }
            })
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.elements[i][j]
    }

    fn take(&mut self, i: usize, j: usize) -> T {
        std::mem::replace(&mut self.elements[i][j], T::zero())
    }

    pub fn transpose(mut self) -> Self {
        Self::from_fn(self.dimension(), |i, j| self.take(j, i))
    }

    pub fn order(&self) -> usize {
        std::iter::repeat(self.clone())
            .scan(SquareMatrix::one(self.dimension()), |matrix_acc, matrix| {
                let matrix_acc_clone = matrix_acc.clone();
                *matrix_acc *= matrix;
                Some(matrix_acc_clone)
            })
            .enumerate()
            .skip(1)
            .find(|(_, matrix)| matrix.is_one())
            .map(|(order, _)| order)
            .unwrap()
    }

    pub fn determinant(&self) -> T
    where
        T: Debug,
    {
        println!("Before permutation");
        let a = T::sum(
            (0..self.dimension())
                .permutations(self.dimension())
                .filter(|permutation| {
                    (0..self.dimension()).all(|i| !self.get(i, permutation[i]).is_zero())
                })
                .map(|permutation| {
                    println!("In permutation: {permutation:?}");
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
                    let product = T::prod(
                        (0..self.dimension())
                            .map(|i| self.get(i, permutation[i]))
                            .cloned(),
                    );
                    if neg_sign {
                        product.neg()
                    } else {
                        product
                    }
                }),
        );
        println!("After permutation");
        a
    }

    // pub fn characteristic_polynomial(&self) -> Polynomial<T> {
    //     SquareMatrix::from_fn(self.dimension(), |i, j| {
    //         alg::AdditiveMagma::sub(
    //             Polynomial::monomial(self.get(i, j).clone(), alg::AdditiveIdentity::zero()),
    //             if i == j {
    //                 Polynomial::monomial(
    //                     alg::MultiplicativeIdentity::one(),
    //                     alg::MultiplicativeIdentity::one(),
    //                 )
    //             } else {
    //                 alg::AdditiveIdentity::zero()
    //             },
    //         )
    //     })
    //     .determinant()
    // }
}

impl<T> std::ops::Neg for SquareMatrix<T>
where
    T: Ring + Clone,
{
    type Output = Self;

    fn neg(mut self) -> Self {
        Self::from_fn(self.dimension(), |i, j| self.take(i, j).neg())
    }
}

impl<T> std::ops::Add for SquareMatrix<T>
where
    T: Ring + Clone,
{
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self {
        assert_eq!(self.dimension(), rhs.dimension());
        Self::from_fn(self.dimension(), |i, j| self.take(i, j).add(rhs.take(i, j)))
    }
}

impl<T> std::ops::Sub for SquareMatrix<T>
where
    T: Ring + Clone,
{
    type Output = Self;

    fn sub(mut self, mut rhs: Self) -> Self {
        assert_eq!(self.dimension(), rhs.dimension());
        Self::from_fn(self.dimension(), |i, j| self.take(i, j).sub(rhs.take(i, j)))
    }
}

impl<T> std::ops::Mul for SquareMatrix<T>
where
    T: Ring + Clone,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.dimension(), rhs.dimension());
        Self::from_fn(self.dimension(), |i, j| {
            T::sum((0..self.dimension()).map(|k| self.get(i, k).clone().mul(rhs.get(k, j).clone())))
        })
    }
}

impl<T> std::ops::AddAssign for SquareMatrix<T>
where
    T: Ring + Clone,
{
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs + rhs;
    }
}

impl<T> std::ops::SubAssign for SquareMatrix<T>
where
    T: Ring + Clone,
{
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs - rhs;
    }
}

impl<T> std::ops::MulAssign for SquareMatrix<T>
where
    T: Ring + Clone,
{
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs * rhs;
    }
}

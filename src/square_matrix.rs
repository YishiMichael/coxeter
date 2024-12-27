use itertools::Itertools;

use super::alg;
use super::polynomial::Polynomial;

#[derive(Clone)]
pub struct SquareMatrix<T>(
    petgraph::matrix_graph::MatrixGraph<(), T, petgraph::Directed, Option<T>, usize>,
);

impl<T> SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    pub fn from_fn<F>(dimension: usize, mut f: F) -> Self
    where
        F: FnMut(usize, usize) -> T,
    {
        Self(petgraph::matrix_graph::MatrixGraph::from_edges(
            (0..dimension)
                .cartesian_product(0..dimension)
                .map(|(i, j)| (i, j, f(i, j))),
        ))
    }

    pub fn zero(dimension: usize) -> Self {
        Self::from_fn(dimension, |_, _| alg::AdditiveIdentity::zero())
    }

    pub fn is_zero(&self) -> bool {
        (0..self.dimension())
            .cartesian_product(0..self.dimension())
            .all(|(i, j)| self.get(i, j).is_zero())
    }

    pub fn one(dimension: usize) -> Self {
        Self::from_fn(dimension, |i, j| {
            if i == j {
                alg::MultiplicativeIdentity::one()
            } else {
                alg::AdditiveIdentity::zero()
            }
        })
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
        self.0.node_count()
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        let indexing = |index| petgraph::visit::NodeIndexable::from_index(&self.0, index);
        self.0.edge_weight(indexing(i), indexing(j))
    }

    pub fn take(&mut self, i: usize, j: usize) -> T {
        let indexing = |index| petgraph::visit::NodeIndexable::from_index(&self.0, index);
        std::mem::replace(
            self.0.edge_weight_mut(indexing(i), indexing(j)),
            alg::AdditiveIdentity::zero(),
        )
    }

    pub fn transpose(mut self) -> Self {
        Self::from_fn(self.dimension(), |i, j| self.take(j, i))
    }

    pub fn order_unipotent(&self) -> usize {
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

    pub fn inverse_one_minus_nilpotent(&self) -> Self {
        std::iter::repeat(self.clone())
            .scan(SquareMatrix::one(self.dimension()), |matrix_acc, matrix| {
                let matrix_acc_clone = matrix_acc.clone();
                *matrix_acc *= matrix;
                Some(matrix_acc_clone)
            })
            .take(self.dimension())
            .fold(SquareMatrix::zero(self.dimension()), std::ops::Add::add)
    }

    pub fn determinant(&self) -> T {
        alg::AdditiveMagma::sum((0..self.dimension()).permutations(self.dimension()).map(
            |permutation| {
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
                let product = alg::MultiplicativeMagma::product(
                    (0..self.dimension())
                        .map(|i| self.get(i, permutation[i]))
                        .cloned(),
                );
                if neg_sign {
                    product.neg()
                } else {
                    product
                }
            },
        ))
    }

    pub fn characteristic_polynomial(&self) -> Polynomial<T> {
        SquareMatrix::from_fn(self.dimension(), |i, j| {
            alg::AdditiveMagma::add(
                if i == j {
                    Polynomial::monomial(
                        alg::MultiplicativeIdentity::one(),
                        alg::MultiplicativeIdentity::one(),
                    )
                } else {
                    alg::AdditiveIdentity::zero()
                },
                alg::AdditiveInverse::neg(Polynomial::monomial(
                    self.get(i, j).clone(),
                    alg::AdditiveIdentity::zero(),
                )),
            )
        })
        .determinant()
    }
}

impl<T> std::ops::Neg for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    type Output = Self;

    fn neg(mut self) -> Self {
        Self::from_fn(self.dimension(), |i, j| self.take(i, j).neg())
    }
}

impl<T> std::ops::Add for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self {
        assert_eq!(self.dimension(), rhs.dimension());
        Self::from_fn(self.dimension(), |i, j| {
            alg::AdditiveMagma::add(self.take(i, j), rhs.take(i, j))
        })
    }
}

impl<T> std::ops::Sub for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    type Output = Self;

    fn sub(mut self, mut rhs: Self) -> Self {
        assert_eq!(self.dimension(), rhs.dimension());
        Self::from_fn(self.dimension(), |i, j| {
            alg::AdditiveMagma::add(self.take(i, j), rhs.take(i, j).neg())
        })
    }
}

impl<T> std::ops::Mul for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.dimension(), rhs.dimension());
        Self::from_fn(self.dimension(), |i, j| {
            alg::AdditiveMagma::sum((0..self.dimension()).map(|k| {
                alg::MultiplicativeMagma::mul(self.get(i, k).clone(), rhs.get(k, j).clone())
            }))
        })
    }
}

impl<T> std::ops::AddAssign for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs + rhs;
    }
}

impl<T> std::ops::SubAssign for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs - rhs;
    }
}

impl<T> std::ops::MulAssign for SquareMatrix<T>
where
    T: alg::Ring + Clone,
{
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs * rhs;
    }
}

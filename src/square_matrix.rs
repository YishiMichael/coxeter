use itertools::Itertools;
// use num::Integer;

#[derive(Clone)]
pub struct SquareMatrix<T> {
    dimension: usize,
    elements: Vec<T>,
}

impl<T> SquareMatrix<T> {
    pub fn from_fn<F>(dimension: usize, mut f: F) -> Self
    where
        F: FnMut((usize, usize)) -> T,
    {
        Self {
            dimension,
            elements: (0..dimension * dimension)
                .map(|index| f(index.div_mod_floor(&dimension)))
                .collect(),
        }
    }

    pub fn map<U, F>(self, f: F) -> SquareMatrix<U>
    where
        F: FnMut(T) -> U,
    {
        SquareMatrix {
            dimension: self.dimension,
            elements: self.elements.into_iter().map(f).collect(),
        }
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn get(&self, (i, j): (usize, usize)) -> &T {
        &self.elements[i * self.dimension + j]
    }
}

impl<T> SquareMatrix<T>
where
    T: Clone,
{
    pub fn embed(dimension: usize, value: T) -> Self {
        Self::from_fn(
            dimension,
            |(i, j)| {
                if i == j {
                    value.clone()
                } else {
                    T::zero()
                }
            },
        )
    }

    fn take(&mut self, (i, j): (usize, usize)) -> T {
        std::mem::replace(&mut self.elements[i * self.dimension + j], T::zero())
    }

    pub fn transpose(mut self) -> Self {
        Self::from_fn(self.dimension, |(i, j)| self.take((j, i)))
    }

    pub fn zero(dimension: usize) -> Self {
        Self::from_fn(dimension, |_| T::zero())
    }

    pub fn is_zero(&self) -> bool {
        self.elements.iter().all(|element| element.is_zero())
    }

    pub fn one(dimension: usize) -> Self {
        Self::from_fn(
            dimension,
            |(i, j)| if i == j { T::one() } else { T::zero() },
        )
    }

    pub fn is_one(&self) -> bool {
        (self.clone() - Self::one(self.dimension)).is_zero()
    }

    pub fn order(&self) -> usize {
        std::iter::repeat(self.clone())
            .scan(SquareMatrix::one(self.dimension), |matrix_acc, matrix| {
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

    pub fn determinant(&self) -> T {
        T::sum(
            (0..self.dimension)
                .permutations(self.dimension)
                .filter(|permutation| {
                    (0..self.dimension).all(|i| !self.get((i, permutation[i])).is_zero())
                })
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
                    let product = T::prod(
                        (0..self.dimension)
                            .map(|i| self.get((i, permutation[i])))
                            .cloned(),
                    );
                    if neg_sign {
                        product.neg()
                    } else {
                        product
                    }
                }),
        )
    }
}

impl<T> std::ops::Neg for SquareMatrix<T>
where
    T: Clone,
{
    type Output = Self;

    fn neg(mut self) -> Self {
        Self::from_fn(self.dimension, |(i, j)| self.take((i, j)).neg())
    }
}

impl<T> std::ops::Add for SquareMatrix<T>
where
    T: Clone,
{
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self {
        assert_eq!(self.dimension, rhs.dimension);
        Self::from_fn(self.dimension, |(i, j)| {
            self.take((i, j)).add(rhs.take((i, j)))
        })
    }
}

impl<T> std::ops::Sub for SquareMatrix<T>
where
    T: Clone,
{
    type Output = Self;

    fn sub(mut self, mut rhs: Self) -> Self {
        assert_eq!(self.dimension, rhs.dimension);
        Self::from_fn(self.dimension, |(i, j)| {
            self.take((i, j)).sub(rhs.take((i, j)))
        })
    }
}

impl<T> std::ops::Mul for SquareMatrix<T>
where
    T: Clone,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.dimension, rhs.dimension);
        Self::from_fn(self.dimension, |(i, j)| {
            T::sum(
                (0..self.dimension).map(|k| self.get((i, k)).clone().mul(rhs.get((k, j)).clone())),
            )
        })
    }
}

impl<T> std::ops::AddAssign for SquareMatrix<T>
where
    T: Clone,
{
    fn add_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs + rhs;
    }
}

impl<T> std::ops::SubAssign for SquareMatrix<T>
where
    T: Clone,
{
    fn sub_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs - rhs;
    }
}

impl<T> std::ops::MulAssign for SquareMatrix<T>
where
    T: Clone,
{
    fn mul_assign(&mut self, rhs: Self) {
        let lhs = std::mem::replace(self, Self::zero(0));
        *self = lhs * rhs;
    }
}

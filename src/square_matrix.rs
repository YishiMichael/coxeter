use feanor_math::{
    matrix::OwnedMatrix,
    ring::{El, RingBase, RingStore},
};
use itertools::Itertools;

#[derive(Clone)]
pub struct SquareMatrixRingBase<R> {
    base_ring: R,
    dimension: usize,
}

impl<R> SquareMatrixRingBase<R> {
    pub fn new(base_ring: R, dimension: usize) -> Self {
        Self {
            base_ring,
            dimension,
        }
    }
}

impl<R> PartialEq for SquareMatrixRingBase<R>
where
    R: RingStore,
{
    fn eq(&self, other: &Self) -> bool {
        self.base_ring.get_ring() == other.base_ring.get_ring() && self.dimension == other.dimension
    }
}

impl<R> RingBase for SquareMatrixRingBase<R>
where
    R: RingStore,
{
    type Element = OwnedMatrix<El<R>>;

    fn clone_el(&self, value: &Self::Element) -> Self::Element {
        value.clone_matrix(&self.base_ring)
    }

    fn add_assign(&self, lhs: &mut Self::Element, mut rhs: Self::Element) {
        (0..self.dimension)
            .cartesian_product(0..self.dimension)
            .for_each(|(i, j)| {
                self.base_ring.add_assign(
                    lhs.at_mut(i, j),
                    std::mem::replace(rhs.at_mut(i, j), self.base_ring.zero()),
                );
            });
    }

    fn negate_inplace(&self, lhs: &mut Self::Element) {
        (0..self.dimension)
            .cartesian_product(0..self.dimension)
            .for_each(|(i, j)| self.base_ring.negate_inplace(lhs.at_mut(i, j)));
    }

    fn mul_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        let els = OwnedMatrix::from_fn(self.dimension, self.dimension, |i, j| {
            self.base_ring.sum(
                (0..self.dimension).map(|k| self.base_ring.mul_ref(lhs.at(i, k), rhs.at(k, j))),
            )
        });
        *lhs = els;
    }

    fn from_int(&self, value: i32) -> Self::Element {
        OwnedMatrix::from_fn(self.dimension, self.dimension, |i, j| {
            if i == j {
                self.base_ring.get_ring().from_int(value)
            } else {
                self.base_ring.zero()
            }
        })
    }

    fn eq_el(&self, lhs: &Self::Element, rhs: &Self::Element) -> bool {
        (0..self.dimension)
            .cartesian_product(0..self.dimension)
            .all(|(i, j)| self.base_ring.eq_el(lhs.at(i, j), rhs.at(i, j)))
    }

    fn is_commutative(&self) -> bool {
        false
    }

    fn is_noetherian(&self) -> bool {
        false
    }

    fn is_approximate(&self) -> bool {
        self.base_ring.get_ring().is_approximate()
    }

    fn dbg_within<'a>(
        &self,
        value: &Self::Element,
        out: &mut std::fmt::Formatter<'a>,
        env: feanor_math::ring::EnvBindingStrength,
    ) -> std::fmt::Result {
        out.write_str(&"[")?;
        let mut row_iter = value.data().row_iter();
        if let Some(first_row) = row_iter.next() {
            out.write_str(&"[")?;
            let mut el_iter = first_row.into_iter();
            if let Some(first_el) = el_iter.next() {
                self.base_ring.get_ring().dbg_within(first_el, out, env)?;
            }
            for el in el_iter {
                out.write_str(&", ")?;
                self.base_ring.get_ring().dbg_within(el, out, env)?;
            }
            out.write_str(&"]")?;
        }
        for row in row_iter {
            out.write_str(&", ")?;
            out.write_str(&"[")?;
            let mut el_iter = row.into_iter();
            if let Some(first_el) = el_iter.next() {
                self.base_ring.get_ring().dbg_within(first_el, out, env)?;
            }
            for el in el_iter {
                out.write_str(&", ")?;
                self.base_ring.get_ring().dbg_within(el, out, env)?;
            }
            out.write_str(&"]")?;
        }
        out.write_str(&"]\n")?;
        Ok(())
    }

    fn characteristic<I>(&self, integer_ring: I) -> Option<El<I>>
    where
        I: RingStore + Copy,
        I::Type: feanor_math::integer::IntegerRing,
    {
        self.base_ring.characteristic(integer_ring)
    }
}

pub fn determinant<R>(
    ring: &impl RingStore<Type = SquareMatrixRingBase<R>>,
    matrix: OwnedMatrix<El<R>>,
) -> El<R>
where
    R: RingStore,
{
    let dimension = ring.get_ring().dimension;
    let base_ring = &ring.get_ring().base_ring;
    base_ring.sum(
        (0..dimension)
            .permutations(dimension)
            .filter(|permutation| {
                (0..dimension).all(|i| !base_ring.is_zero(matrix.at(i, permutation[i])))
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
                let product = base_ring
                    .prod((0..dimension).map(|i| base_ring.clone_el(matrix.at(i, permutation[i]))));
                if neg_sign {
                    base_ring.negate(product)
                } else {
                    product
                }
            }),
    )
}

use std::collections::HashMap;

pub trait Order: num::Integer + num::Unsigned + std::hash::Hash + Copy {}

impl Order for u32 {}

pub trait RingElement:
    num::Zero
    + std::ops::Add
    + std::ops::AddAssign
    + std::ops::Sub
    + std::ops::SubAssign
    + std::ops::Neg
    + num::One
    + std::ops::Mul
    + std::ops::MulAssign
{
}

impl RingElement for i32 {}

// Ring structure R[zeta_n]
// Basis of coefficients: see Thomas, Integral Bases for Subfields of Cyclotomic Fields
// Sparse representation. Entries in `coeff_map` must be nonzero.
#[derive(Clone)]
pub struct CyclotomicNumber<C, D = u32>
where
    D: Order,
{
    order: D,
    coeff_map: HashMap<D, C>,
}

impl<C, D> CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
    pub fn root_of_unity(order: D, exponent: D) -> Self {
        todo!()
    }
}

impl<C, D> num::Zero for CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
    fn zero() -> Self {
        Self {
            order: D::one(),
            coeff_map: HashMap::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.coeff_map.is_empty()
    }
}

impl<C, D> num::One for CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
    fn one() -> Self {
        Self {
            order: D::one(),
            coeff_map: {
                let mut coeff_map = HashMap::new();
                coeff_map.insert(D::one(), C::one());
                coeff_map
            },
        }
    }
}

impl<C, D> std::ops::AddAssign for CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
    fn add_assign(&mut self, rhs: Self) {
        let order = self.order.lcm(&rhs.order);
        let scalar = order / self.order;
        let rhs_scalar = order / rhs.order;
        self.order = order;
        self.coeff_map = std::mem::take(&mut self.coeff_map)
            .into_iter()
            .map(|(exponent, coeff)| (scalar * exponent, coeff))
            .collect();
        rhs.coeff_map.into_iter().for_each(|(exponent, coeff)| {
            *self
                .coeff_map
                .entry(rhs_scalar * exponent)
                .or_insert(C::zero()) += coeff
        });
    }
}

impl<C, D> std::ops::Add for CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<C, D> std::ops::Mul for CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl<C, D> RingElement for CyclotomicNumber<C, D>
where
    C: RingElement,
    D: Order,
{
}

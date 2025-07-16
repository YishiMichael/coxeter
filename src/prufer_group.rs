use super::alg::*;

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
pub struct PruferGroup<U> {
    numer: U,
    denom: U,
}

impl<U> PruferGroup<U>
where
    U: num::Integer + num::Unsigned,
{
    #[inline]
    pub fn new(numer: U, denom: U) -> Self {
        assert!(!denom.is_zero());
        Self { numer, denom }
    }

    #[inline]
    pub fn numer(&self) -> &U {
        &self.numer
    }

    #[inline]
    pub fn denom(&self) -> &U {
        &self.denom
    }
}

impl<U> AdditiveMagma for PruferGroup<U>
where
    U: num::Integer + num::Unsigned,
{
    #[inline]
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |lhs, rhs| {
            let denom = lhs.denom.lcm(&rhs.denom);
            let numer = (lhs.numer * denom.div_floor(&lhs.denom)
                + rhs.numer * denom.div_floor(&rhs.denom))
            .mod_floor(&denom);
            Self { numer, denom }
        })
    }
}

impl<U> AdditiveIdentity for PruferGroup<U>
where
    U: num::Integer + num::Unsigned,
{
    #[inline]
    fn zero() -> Self {
        Self::new(U::zero(), U::one())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.numer.is_zero()
    }
}

impl<U> AdditiveInverse for PruferGroup<U>
where
    U: Clone + num::Integer + num::Unsigned,
{
    #[inline]
    fn neg(self) -> Self {
        Self::new(self.denom.clone() - self.numer, self.denom)
    }
}

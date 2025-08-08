use feanor_math::{
    homomorphism::CanHomFrom,
    ring::{El, RingBase, RingStore},
};
use itertools::Itertools;

#[derive(Clone)]
pub struct CyclotomicRingBase<R> {
    base_ring: R,
}

impl<R> CyclotomicRingBase<R> {
    pub fn new(base_ring: R) -> Self {
        Self { base_ring }
    }
}

impl<R> PartialEq for CyclotomicRingBase<R>
where
    R: RingStore,
{
    fn eq(&self, other: &Self) -> bool {
        self.base_ring.get_ring() == other.base_ring.get_ring()
    }
}

impl<R> RingBase for CyclotomicRingBase<R>
where
    R: RingStore,
{
    type Element = Vec<El<R>>;

    fn clone_el(&self, value: &Self::Element) -> Self::Element {
        value
            .iter()
            .map(|el| self.base_ring.get_ring().clone_el(el))
            .collect()
    }

    fn add_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        let lhs_order = lhs.len();
        let rhs_order = rhs.len();
        let order = lcm(lhs_order, rhs_order);
        let mut els = std::iter::repeat_with(|| self.base_ring.zero())
            .take(order)
            .collect_vec();
        els.iter_mut()
            .step_by(order / lhs_order)
            .zip_eq(lhs.drain(..))
            .for_each(|(el, lhs_el)| self.base_ring.add_assign(el, lhs_el));
        els.iter_mut()
            .step_by(order / rhs_order)
            .zip_eq(rhs.into_iter())
            .for_each(|(el, rhs_el)| self.base_ring.add_assign(el, rhs_el));
        *lhs = els;
    }

    fn negate_inplace(&self, lhs: &mut Self::Element) {
        lhs.iter_mut()
            .for_each(|el| self.base_ring.negate_inplace(el));
    }

    fn mul_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        let lhs_order = lhs.len();
        let rhs_order = rhs.len();
        let order = lcm(lhs_order, rhs_order);
        let lhs_step = order / lhs_order;
        let rhs_step = order / rhs_order;
        let mut els = std::iter::repeat_with(|| self.base_ring.zero())
            .take(order)
            .collect_vec();
        lhs.iter()
            .enumerate()
            .filter(|(_, lhs_el)| !self.base_ring.is_zero(lhs_el))
            .cartesian_product(
                rhs.iter()
                    .enumerate()
                    .filter(|(_, rhs_el)| !self.base_ring.is_zero(rhs_el)),
            )
            .for_each(|((lhs_index, lhs_el), (rhs_index, rhs_el))| {
                self.base_ring.add_assign(
                    &mut els[(lhs_index * lhs_step + rhs_index * rhs_step) % order],
                    self.base_ring.mul_ref(lhs_el, rhs_el),
                );
            });
        *lhs = els;
    }

    fn from_int(&self, value: i32) -> Self::Element {
        Vec::from([self.base_ring.get_ring().from_int(value)])
    }

    fn is_zero(&self, value: &Self::Element) -> bool {
        let order = value.len() as u64;
        let balancer = self.prod(
            (2..)
                .scan(order, |n, p| {
                    (*n != 1).then(|| {
                        std::iter::from_fn(|| (*n % p == 0).then(|| *n /= p))
                            .last()
                            .map(|()| p)
                    })
                })
                .filter_map(|option_p| option_p)
                .map(|p| {
                    let mut els = std::iter::repeat_with(|| self.base_ring.one())
                        .take(p as usize)
                        .collect_vec();
                    self.base_ring
                        .sub_assign(&mut els[0], self.base_ring.get_ring().from_int(p as i32));
                    els
                }),
        );
        self.mul_ref_snd(balancer, value)
            .iter()
            .all(|el| self.base_ring.is_zero(el))
    }

    fn eq_el(&self, lhs: &Self::Element, rhs: &Self::Element) -> bool {
        self.is_zero(&self.sub_ref(lhs, rhs))
    }

    fn is_commutative(&self) -> bool {
        self.base_ring.is_commutative()
    }

    fn is_noetherian(&self) -> bool {
        false
    }

    fn is_approximate(&self) -> bool {
        false
    }

    fn dbg_within<'a>(
        &self,
        value: &Self::Element,
        out: &mut std::fmt::Formatter<'a>,
        env: feanor_math::ring::EnvBindingStrength,
    ) -> std::fmt::Result {
        let poly_ring = feanor_math::rings::poly::dense_poly::DensePolyRing::new(
            feanor_math::ring::RingRef::new(self.base_ring.get_ring()),
            "X",
        );
        let el = feanor_math::rings::poly::PolyRingStore::from_terms(
            &poly_ring,
            value
                .iter()
                .enumerate()
                .map(|(index, el)| (self.base_ring.clone_el(el), index)),
        );
        feanor_math::rings::poly::generic_impls::dbg_poly(
            poly_ring.get_ring(),
            &el,
            out,
            format!("zeta_{}", value.len()).as_str(),
            env,
        )
    }

    fn characteristic<I>(&self, integer_ring: I) -> Option<El<I>>
    where
        I: RingStore + Copy,
        I::Type: feanor_math::integer::IntegerRing,
    {
        self.base_ring.characteristic(integer_ring)
    }
}

fn lcm(a: usize, b: usize) -> usize {
    debug_assert!(a != 0);
    debug_assert!(b != 0);
    let mut a0 = a;
    let mut b0 = b;
    while b0 != 0 {
        (a0, b0) = (b0, a0 % b0);
    }
    a / a0 * b
}

pub fn cyclotomic_root_of_unity<R>(
    ring: &impl RingStore<Type = CyclotomicRingBase<R>>,
    order: u64,
    exponent: u64,
) -> <CyclotomicRingBase<R> as RingBase>::Element
where
    R: RingStore,
{
    debug_assert!(order != 0);
    let base_ring = &ring.get_ring().base_ring;
    let mut els = std::iter::repeat_with(|| base_ring.zero())
        .take(order as usize)
        .collect_vec();
    base_ring.add_assign(
        &mut els[exponent as usize % order as usize],
        base_ring.one(),
    );
    els
}

pub fn cyclotomic_numeric_embed_into<R>(
    ring: &impl RingStore<Type = CyclotomicRingBase<R>>,
    value: <CyclotomicRingBase<R> as RingBase>::Element,
) -> feanor_math::rings::float_complex::Complex64El
where
    R: RingStore,
    feanor_math::rings::approx_real::float::Real64Base: CanHomFrom<R::Type>,
{
    const RR: feanor_math::rings::approx_real::float::Real64 =
        feanor_math::rings::approx_real::float::Real64::RING;
    const CC: feanor_math::rings::float_complex::Complex64 =
        feanor_math::rings::float_complex::Complex64::RING;
    let order = value.len();
    CC.sum(value.into_iter().enumerate().map(|(index, el)| {
        CC.mul(
            CC.from_f64(
                RR.get_ring().map_in(
                    ring.get_ring().base_ring.get_ring(),
                    el,
                    &RR.get_ring()
                        .has_canonical_hom(ring.get_ring().base_ring.get_ring())
                        .unwrap(),
                ),
            ),
            CC.root_of_unity(index as i64, order as i64),
        )
    }))
}

#[cfg(test)]
mod test {
    use super::{cyclotomic_root_of_unity, CyclotomicRingBase};
    use feanor_math::ring::{RingStore, RingValue};

    #[test]
    fn test_cyclotomic() {
        let cr = RingValue::from(CyclotomicRingBase::new(
            feanor_math::primitive_int::StaticRing::<i64>::default(),
        ));
        assert!(cr.is_zero(&cr.sum([
            cyclotomic_root_of_unity(&cr, 3, 0),
            cr.negate(cyclotomic_root_of_unity(&cr, 2, 0)),
        ])));
        assert!(cr.is_zero(&cr.sum([
            cyclotomic_root_of_unity(&cr, 3, 1),
            cyclotomic_root_of_unity(&cr, 3, 2),
            cyclotomic_root_of_unity(&cr, 3, 0),
        ])));
        assert!(cr.is_zero(&cr.sum([
            cyclotomic_root_of_unity(&cr, 3, 1),
            cyclotomic_root_of_unity(&cr, 9, 6),
            cyclotomic_root_of_unity(&cr, 3, 0),
        ])));
        assert!(cr.is_zero(&cr.sum([
            cyclotomic_root_of_unity(&cr, 2, 0),
            cyclotomic_root_of_unity(&cr, 6, 2),
            cyclotomic_root_of_unity(&cr, 3, 2),
        ])));
        assert!(cr.is_zero(&cr.sum([
            cyclotomic_root_of_unity(&cr, 5, 1),
            cyclotomic_root_of_unity(&cr, 5, 2),
            cyclotomic_root_of_unity(&cr, 5, 3),
            cyclotomic_root_of_unity(&cr, 5, 4),
            cr.negate(cyclotomic_root_of_unity(&cr, 3, 1)),
            cr.negate(cyclotomic_root_of_unity(&cr, 3, 2))
        ])));
    }
}

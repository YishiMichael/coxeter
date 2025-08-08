use feanor_math::{
    matrix::OwnedMatrix,
    ring::{El, RingBase, RingStore, RingValue},
    rings::poly::PolyRingStore,
};
use itertools::Itertools;

use super::dynkin_diagram::DynkinDiagramType;
use super::square_matrix::SquareMatrixRingBase;

pub enum Label {
    Implicit,
    Explicit(u64),
}

#[derive(Debug, PartialEq)]
pub struct CoxeterGroupInfo {
    pub coxeter_number: u64,
    pub order: u64,
    pub exponents: Vec<u64>,
    pub degrees: Vec<u64>,
}

pub struct CoxeterDiagram<N>(petgraph::graph::UnGraph<N, Label>);

impl<N> CoxeterDiagram<N> {
    pub fn new<I>(rank: usize, labels: I) -> Self
    where
        N: Default,
        I: IntoIterator<Item = ((usize, usize), Label)>,
    {
        let mut graph = petgraph::graph::UnGraph::new_undirected();
        let nodes = (0..rank)
            .map(|_| graph.add_node(N::default()))
            .collect::<Vec<_>>();
        labels.into_iter().for_each(|((i, j), label)| {
            graph.add_edge(nodes[i], nodes[j], label);
        });
        Self(graph)
    }

    pub fn rank(&self) -> usize {
        self.0.node_count()
    }

    fn cyclotomic_ring(
        &self,
    ) -> RingValue<CyclotomicRingBase<feanor_math::primitive_int::StaticRing<i64>>> {
        RingValue::from(CyclotomicRingBase::new(
            feanor_math::primitive_int::StaticRing::<i64>::default(),
        ))
    }

    fn matrix_ring(
        &self,
    ) -> RingValue<
        SquareMatrixRingBase<
            RingValue<CyclotomicRingBase<feanor_math::primitive_int::StaticRing<i64>>>,
        >,
    > {
        RingValue::from(SquareMatrixRingBase::new(
            self.cyclotomic_ring(),
            self.rank(),
        ))
    }

    fn get_edge(&self, (i, j): (usize, usize)) -> Option<&Label> {
        self.0
            .find_edge(
                petgraph::visit::NodeIndexable::from_index(&self.0, i),
                petgraph::visit::NodeIndexable::from_index(&self.0, j),
            )
            .map(|edge_index| self.0.edge_weight(edge_index).unwrap())
    }

    // pub fn coxeter_matrix(&self) -> OwnedMatrix<u64> {
    //     let rank = self.rank();
    //     OwnedMatrix::from_fn(rank, rank, |i, j| {
    //         self.get_edge((i, j))
    //             .map(|label| match label {
    //                 Label::Implicit => 3,
    //                 Label::Explicit(label) => *label,
    //             })
    //             .unwrap_or_else(|| if i == j { 1 } else { 2 })
    //     })
    // }

    pub fn schlafli_matrix(
        &self,
    ) -> <SquareMatrixRingBase<
        RingValue<CyclotomicRingBase<feanor_math::primitive_int::StaticRing<i64>>>,
    > as RingBase>::Element {
        let rank = self.rank();
        let cr = self.cyclotomic_ring();
        OwnedMatrix::from_fn(rank, rank, |i, j| {
            self.get_edge((i, j))
                .map(|label| match label {
                    Label::Implicit => cr.get_ring().from_int(-1),
                    Label::Explicit(label) => cr.negate(cr.add(
                        root_of_unity(&cr, 2 * label, 1),
                        root_of_unity(&cr, 2 * label, 2 * label - 1),
                    )),
                })
                .unwrap_or_else(|| {
                    if i == j {
                        cr.get_ring().from_int(2)
                    } else {
                        cr.zero()
                    }
                })
        })
    }

    pub fn coxeter_element(
        &self,
    ) -> <SquareMatrixRingBase<
        RingValue<CyclotomicRingBase<feanor_math::primitive_int::StaticRing<i64>>>,
    > as RingBase>::Element {
        let rank = self.rank();
        let cr = self.cyclotomic_ring();
        let mr = self.matrix_ring();
        let mut schlafli_matrix = self.schlafli_matrix();
        let strict_upper_triangular_matrix = OwnedMatrix::from_fn(rank, rank, |i, j| {
            if i < j {
                std::mem::replace(schlafli_matrix.at_mut(i, j), cr.zero())
            } else {
                cr.zero()
            }
        });
        let strict_lower_triangular_matrix = OwnedMatrix::from_fn(rank, rank, |i, j| {
            if i > j {
                std::mem::replace(schlafli_matrix.at_mut(i, j), cr.zero())
            } else {
                cr.zero()
            }
        });
        mr.mul(
            mr.negate(mr.add(mr.one(), strict_lower_triangular_matrix)),
            mr.sum(
                std::iter::successors(Some(mr.one()), |matrix| {
                    Some(mr.negate(mr.mul_ref(&strict_upper_triangular_matrix, matrix)))
                })
                .take(rank),
            ),
        )
    }

    pub fn coxeter_group_info(&self) -> CoxeterGroupInfo {
        let rank = self.rank();
        let cr = self.cyclotomic_ring();
        let mr = self.matrix_ring();
        let pr = feanor_math::rings::poly::dense_poly::DensePolyRing::new(cr.clone(), "lambda");
        let mut coxeter_element = self.coxeter_element();
        let coxeter_number = std::iter::successors(Some(mr.one()), |matrix| {
            Some(mr.mul_ref(&coxeter_element, matrix))
        })
        .enumerate()
        .skip(1)
        .find(|(_, matrix)| mr.is_one(matrix))
        .map(|(order, _)| order)
        .unwrap() as u64;
        let polynomial_matrix = OwnedMatrix::from_fn(rank, rank, |i, j| {
            if i == j {
                pr.from_terms([
                    (
                        cr.negate(std::mem::replace(coxeter_element.at_mut(i, j), cr.zero())),
                        0,
                    ),
                    (cr.one(), 1),
                ])
            } else {
                pr.from_terms([(
                    cr.negate(std::mem::replace(coxeter_element.at_mut(i, j), cr.zero())),
                    0,
                )])
            }
        });
        let characteristic_polynomial = pr.sum(
            (0..rank)
                .permutations(rank)
                .filter(|permutation| {
                    (0..rank).all(|i| !pr.is_zero(polynomial_matrix.at(i, permutation[i])))
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
                    let product = pr.prod(
                        (0..rank).map(|i| pr.clone_el(polynomial_matrix.at(i, permutation[i]))),
                    );
                    if neg_sign {
                        pr.negate(product)
                    } else {
                        product
                    }
                }),
        );
        let exponents = (0..coxeter_number)
            .scan(characteristic_polynomial, |polynomial, exponent| {
                (!pr.is_one(&polynomial)).then(|| {
                    let monomial = pr.from_terms([
                        (cr.negate(root_of_unity(&cr, coxeter_number, exponent)), 0),
                        (cr.one(), 1),
                    ]);
                    std::iter::from_fn(|| {
                        let (quo, rem) = pr.div_rem_monic(pr.clone_el(polynomial), &monomial);
                        if pr.is_zero(&rem) {
                            *polynomial = quo;
                            Some(exponent)
                        } else {
                            None
                        }
                    })
                    .collect_vec()
                })
            })
            .flatten()
            .collect_vec();
        let degrees = exponents.iter().map(|&exponent| exponent + 1).collect_vec();
        CoxeterGroupInfo {
            coxeter_number,
            order: degrees.iter().product(),
            exponents,
            degrees,
        }
    }
}

pub enum CoxeterDiagramType {
    A(usize),
    B(usize),
    C(usize),
    D(usize),
    E(usize),
    F(usize),
    G(usize),
    H(usize),
    I2(u64),
}

impl From<CoxeterDiagramType> for CoxeterDiagram<()> {
    fn from(value: CoxeterDiagramType) -> Self {
        match value {
            CoxeterDiagramType::A(rank @ 1..) => Self::new(
                rank,
                (0..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit)),
            ),
            CoxeterDiagramType::B(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([((rank - 2, rank - 1), Label::Explicit(4))]),
            ),
            CoxeterDiagramType::C(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([((rank - 2, rank - 1), Label::Explicit(4))]),
            ),
            CoxeterDiagramType::D(rank @ 3..) => Self::new(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([
                        ((rank - 3, rank - 2), Label::Implicit),
                        ((rank - 3, rank - 1), Label::Implicit),
                    ]),
            ),
            CoxeterDiagramType::E(rank @ 4..=8) => Self::new(
                rank,
                (2..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([((0, 2), Label::Implicit), ((1, 3), Label::Implicit)]),
            ),
            CoxeterDiagramType::F(rank @ 3..=4) => Self::new(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([
                        ((rank - 3, rank - 2), Label::Explicit(4)),
                        ((rank - 2, rank - 1), Label::Implicit),
                    ]),
            ),
            CoxeterDiagramType::G(rank @ 2..=2) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([((rank - 2, rank - 1), Label::Explicit(6))]),
            ),
            CoxeterDiagramType::H(rank @ 2..=4) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Label::Implicit))
                    .chain([((rank - 2, rank - 1), Label::Explicit(5))]),
            ),
            CoxeterDiagramType::I2(label @ 3..) => Self::new(2, [((0, 1), Label::Explicit(label))]),
            _ => panic!(),
        }
    }
}

impl From<DynkinDiagramType> for CoxeterDiagramType {
    fn from(value: DynkinDiagramType) -> Self {
        match value {
            DynkinDiagramType::A(rank) => Self::A(rank),
            DynkinDiagramType::B(rank) => Self::B(rank),
            DynkinDiagramType::C(rank) => Self::C(rank),
            DynkinDiagramType::D(rank) => Self::D(rank),
            DynkinDiagramType::E(rank) => Self::E(rank),
            DynkinDiagramType::F(rank) => Self::F(rank),
            DynkinDiagramType::G(rank) => Self::G(rank),
        }
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

#[derive(Clone)]
pub struct CyclotomicRingBase<R> {
    base_ring: R,
}

impl<R> CyclotomicRingBase<R> {
    fn new(base_ring: R) -> Self {
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

fn root_of_unity<R>(
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

#[cfg(test)]
mod test {
    use feanor_math::ring::{RingStore, RingValue};

    use crate::coxeter_diagram::CoxeterGroupInfo;

    use super::root_of_unity;
    use super::CoxeterDiagram;
    use super::CoxeterDiagramType;
    use super::CyclotomicRingBase;

    #[test]
    fn test_cyclotomic() {
        let cr = RingValue::from(CyclotomicRingBase::new(
            feanor_math::primitive_int::StaticRing::<i64>::default(),
        ));
        assert!(cr.is_zero(&cr.sum([
            root_of_unity(&cr, 3, 0),
            cr.negate(root_of_unity(&cr, 2, 0)),
        ])));
        assert!(cr.is_zero(&cr.sum([
            root_of_unity(&cr, 3, 1),
            root_of_unity(&cr, 3, 2),
            root_of_unity(&cr, 3, 0),
        ])));
        assert!(cr.is_zero(&cr.sum([
            root_of_unity(&cr, 3, 1),
            root_of_unity(&cr, 9, 6),
            root_of_unity(&cr, 3, 0),
        ])));
        assert!(cr.is_zero(&cr.sum([
            root_of_unity(&cr, 2, 0),
            root_of_unity(&cr, 6, 2),
            root_of_unity(&cr, 3, 2),
        ])));
        assert!(cr.is_zero(&cr.sum([
            root_of_unity(&cr, 5, 1),
            root_of_unity(&cr, 5, 2),
            root_of_unity(&cr, 5, 3),
            root_of_unity(&cr, 5, 4),
            cr.negate(root_of_unity(&cr, 3, 1)),
            cr.negate(root_of_unity(&cr, 3, 2))
        ])));
    }

    #[test]
    fn test_a3() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::A(3)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 4,
                order: 24,
                exponents: vec![1, 2, 3],
                degrees: vec![2, 3, 4],
            },
        );
    }

    #[test]
    fn test_b3() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::B(3)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 6,
                order: 48,
                exponents: vec![1, 3, 5],
                degrees: vec![2, 4, 6],
            },
        );
    }

    #[test]
    fn test_d4() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::D(4)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 6,
                order: 192,
                exponents: vec![1, 3, 3, 5],
                degrees: vec![2, 4, 4, 6],
            },
        );
    }

    #[test]
    fn test_h3() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::H(3)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 10,
                order: 120,
                exponents: vec![1, 5, 9],
                degrees: vec![2, 6, 10],
            },
        );
    }

    #[test]
    fn test_h4() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::H(4)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 30,
                order: 14400,
                exponents: vec![1, 11, 19, 29],
                degrees: vec![2, 12, 20, 30],
            },
        );
    }

    #[test]
    fn test_g2() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::G(2)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 6,
                order: 12,
                exponents: vec![1, 5],
                degrees: vec![2, 6],
            },
        );
    }

    #[test]
    fn test_f4() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::F(4)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 12,
                order: 1152,
                exponents: vec![1, 5, 7, 11],
                degrees: vec![2, 6, 8, 12],
            },
        );
    }

    #[test]
    fn test_e6() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(6)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 12,
                order: 51840,
                exponents: vec![1, 4, 5, 7, 8, 11],
                degrees: vec![2, 5, 6, 8, 9, 12],
            },
        );
    }

    #[test]
    fn test_e7() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(7)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 18,
                order: 2903040,
                exponents: vec![1, 5, 7, 9, 11, 13, 17],
                degrees: vec![2, 6, 8, 10, 12, 14, 18],
            },
        );
    }

    #[test]
    fn test_e8() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(8)).coxeter_group_info(),
            CoxeterGroupInfo {
                coxeter_number: 30,
                order: 696729600,
                exponents: vec![1, 7, 11, 13, 17, 19, 23, 29],
                degrees: vec![2, 8, 12, 14, 18, 20, 24, 30],
            },
        );
    }
}

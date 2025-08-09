use super::cyclotomic::{
    cyclotomic_numeric_embed_into, cyclotomic_root_of_unity, CyclotomicRingBase,
};
use super::dynkin_diagram::DynkinDiagramType;
use super::square_matrix::{determinant, SquareMatrixRingBase};
use feanor_math::{
    matrix::OwnedMatrix,
    ring::{RingBase, RingStore, RingValue},
    rings::poly::PolyRingStore,
};
use itertools::Itertools;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoxeterGroupType {
    Elliptic,
    Parabolic,
    Hyperbolic,
}

#[derive(Debug, PartialEq)]
pub struct CoxeterGroupInfo {
    pub coxeter_number: u64,
    pub order: u64,
    pub exponents: Vec<u64>,
    pub degrees: Vec<u64>,
}

#[derive(Clone, Debug)]
pub struct CoxeterDiagram(petgraph::graph::UnGraph<(), u64>);

impl CoxeterDiagram {
    pub fn new<I>(rank: usize, labels: I) -> Self
    where
        I: IntoIterator<Item = ((usize, usize), u64)>,
    {
        let mut graph = petgraph::graph::UnGraph::new_undirected();
        (0..rank).for_each(|_| {
            graph.add_node(());
        });
        labels.into_iter().for_each(|((i, j), label)| {
            graph.add_edge(
                petgraph::graph::NodeIndex::new(i),
                petgraph::graph::NodeIndex::new(j),
                label,
            );
        });
        Self(graph)
    }

    fn get_edge(&self, i: usize, j: usize) -> Option<&u64> {
        self.0
            .find_edge(
                petgraph::graph::NodeIndex::new(i),
                petgraph::graph::NodeIndex::new(j),
            )
            .map(|edge_index| self.0.edge_weight(edge_index).unwrap())
    }

    fn add_node(&mut self) {
        self.0.add_node(());
    }

    fn set_edge(&mut self, i: usize, j: usize, label: u64) {
        self.0.update_edge(
            petgraph::graph::NodeIndex::new(i),
            petgraph::graph::NodeIndex::new(j),
            label,
        );
    }

    fn is_isomorphic(&self, other: &Self) -> bool {
        petgraph::algo::is_isomorphic_matching(&self.0, &other.0, PartialEq::eq, PartialEq::eq)
    }

    fn is_isomorphic_subgraph(&self, other: &Self) -> bool {
        petgraph::algo::is_isomorphic_subgraph_matching(
            &self.0,
            &other.0,
            PartialEq::eq,
            PartialEq::eq,
        )
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

    pub fn schlafli_matrix(
        &self,
    ) -> <SquareMatrixRingBase<
        RingValue<CyclotomicRingBase<feanor_math::primitive_int::StaticRing<i64>>>,
    > as RingBase>::Element {
        let rank = self.rank();
        let cr = self.cyclotomic_ring();
        OwnedMatrix::from_fn(rank, rank, |i, j| {
            self.get_edge(i, j)
                .map(|&label| {
                    if label == 3 {
                        cr.neg_one()
                    } else {
                        cr.negate(cr.add(
                            cyclotomic_root_of_unity(&cr, 2 * label, 1),
                            cyclotomic_root_of_unity(&cr, 2 * label, 2 * label - 1),
                        ))
                    }
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

    pub fn coxeter_group_type(&self) -> CoxeterGroupType {
        const CC: feanor_math::rings::float_complex::Complex64 =
            feanor_math::rings::float_complex::Complex64::RING;
        let rank = self.rank();
        let cr = self.cyclotomic_ring();
        let schlafli_matrix = self.schlafli_matrix();
        let mut coxeter_group_type = CoxeterGroupType::Elliptic;
        (1..=rank).for_each(|n| {
            let numeric = CC.re(cyclotomic_numeric_embed_into(
                &cr,
                determinant(
                    &RingValue::from(SquareMatrixRingBase::new(cr.clone(), n)),
                    OwnedMatrix::from_fn(n, n, |i, j| cr.clone_el(schlafli_matrix.at(i, j))),
                ),
            ));
            if numeric > 16.0 * std::f64::EPSILON {
            } else if numeric >= -16.0 * std::f64::EPSILON {
                if coxeter_group_type == CoxeterGroupType::Elliptic {
                    coxeter_group_type = CoxeterGroupType::Parabolic;
                }
            } else {
                if coxeter_group_type == CoxeterGroupType::Elliptic
                    || coxeter_group_type == CoxeterGroupType::Parabolic
                {
                    coxeter_group_type = CoxeterGroupType::Hyperbolic;
                }
            }
        });
        coxeter_group_type
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
        let characteristic_polynomial = determinant(
            &RingValue::from(SquareMatrixRingBase::new(pr.clone(), rank)),
            OwnedMatrix::from_fn(rank, rank, |i, j| {
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
            }),
        );
        let exponents = (0..coxeter_number)
            .scan(characteristic_polynomial, |polynomial, exponent| {
                (!pr.is_one(&polynomial)).then(|| {
                    let monomial = pr.from_terms([
                        (
                            cr.negate(cyclotomic_root_of_unity(&cr, coxeter_number, exponent)),
                            0,
                        ),
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

    pub fn enumerate_connected_acyclic(
        rank: usize,
        max_depth: usize,
    ) -> Vec<(Self, CoxeterGroupType)> {
        let mut coxeter_diagrams = Vec::new();
        if rank == 1 {
            coxeter_diagrams.push((Self::new(1, std::iter::empty()), CoxeterGroupType::Elliptic));
        } else {
            let mut elliptic_coxeter_diagrams = Vec::new();
            let mut parabolic_coxeter_diagrams = Vec::new();
            let mut hyperbolic_coxeter_diagrams = Vec::new();
            for (coxeter_diagram, coxeter_group_type) in
                Self::enumerate_connected_acyclic(rank - 1, max_depth)
            {
                match coxeter_group_type {
                    CoxeterGroupType::Elliptic => elliptic_coxeter_diagrams.push(coxeter_diagram),
                    CoxeterGroupType::Parabolic => parabolic_coxeter_diagrams.push(coxeter_diagram),
                    CoxeterGroupType::Hyperbolic => {
                        hyperbolic_coxeter_diagrams.push(coxeter_diagram)
                    }
                }
            }

            let mut visited_coxeter_diagrams = Vec::new();
            'out: for elliptic in &elliptic_coxeter_diagrams {
                for i in 0..rank - 1 {
                    for label in 3.. {
                        let mut coxeter_diagram = elliptic.clone();
                        coxeter_diagram.add_node();
                        coxeter_diagram.set_edge(i, rank - 1, label);

                        if visited_coxeter_diagrams
                            .iter()
                            .any(|visited| coxeter_diagram.is_isomorphic(visited))
                        {
                            continue;
                        }
                        visited_coxeter_diagrams.push(coxeter_diagram.clone());

                        if hyperbolic_coxeter_diagrams
                            .iter()
                            .any(|hyperbolic| hyperbolic.is_isomorphic_subgraph(&coxeter_diagram))
                        {
                            break;
                        } else if parabolic_coxeter_diagrams
                            .iter()
                            .any(|parabolic| parabolic.is_isomorphic_subgraph(&coxeter_diagram))
                        {
                            coxeter_diagrams.push((coxeter_diagram, CoxeterGroupType::Hyperbolic));
                            break;
                        } else {
                            let coxeter_group_type = coxeter_diagram.coxeter_group_type();
                            coxeter_diagrams.push((coxeter_diagram, coxeter_group_type));
                            if coxeter_diagrams.len() >= max_depth {
                                break 'out;
                            }
                            if coxeter_group_type == CoxeterGroupType::Hyperbolic {
                                break;
                            }
                            continue;
                        }
                    }
                }
            }
        }
        dbg!(&coxeter_diagrams);
        coxeter_diagrams
        // let mut triplet = Triplet {
        //     elliptic: Vec::new(),
        //     parabolic: Vec::new(),
        //     hyperbolic: Vec::new(),
        // };
        // if n == 1 {
        //     triplet.elliptic.push(SchlafliSymbol::empty()); // {}
        // } else {
        //     let prev_triplet = Self::enumerate_triplet(n - 1);
        //     for elliptic in &prev_triplet.elliptic {
        //         for label in 3.. {
        //             let schlafli_symbol = elliptic.add_label(label);
        //             if prev_triplet.hyperbolic.iter().any(|hyperbolic| {
        //                 hyperbolic
        //                     .coxeter_graph()
        //                     .is_isomorphic_subgraph(schlafli_symbol.coxeter_graph())
        //             }) {
        //                 break;
        //             } else if prev_triplet.parabolic.iter().any(|parabolic| {
        //                 parabolic
        //                     .coxeter_graph()
        //                     .is_isomorphic_subgraph(schlafli_symbol.coxeter_graph())
        //             }) {
        //                 triplet.hyperbolic.push(schlafli_symbol);
        //                 break;
        //             } else if schlafli_symbol.coxeter_graph().is_positive_definite() {
        //                 triplet.elliptic.push(schlafli_symbol);
        //                 continue;
        //             } else if schlafli_symbol.coxeter_graph().is_positive_semidefinite() {
        //                 triplet.parabolic.push(schlafli_symbol);
        //                 continue;
        //             } else {
        //                 triplet.hyperbolic.push(schlafli_symbol);
        //                 break;
        //             }
        //         }
        //     }
        // }
        // triplet
    }

    pub fn enumerate_finite(rank: usize, max_depth: usize) -> Vec<Self> {
        Self::enumerate_connected_acyclic(rank, max_depth)
            .into_iter()
            .filter_map(|(coxeter_diagram, coxeter_group_type)| {
                (coxeter_group_type == CoxeterGroupType::Elliptic).then_some(coxeter_diagram)
            })
            .collect()
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

impl From<CoxeterDiagramType> for CoxeterDiagram {
    fn from(value: CoxeterDiagramType) -> Self {
        match value {
            CoxeterDiagramType::A(rank @ 1..) => {
                Self::new(rank, (0..rank).tuple_windows().map(|(i, j)| ((i, j), 3)))
            }
            CoxeterDiagramType::B(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((rank - 2, rank - 1), 4)]),
            ),
            CoxeterDiagramType::C(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((rank - 2, rank - 1), 4)]),
            ),
            CoxeterDiagramType::D(rank @ 3..) => Self::new(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((rank - 3, rank - 2), 3), ((rank - 3, rank - 1), 3)]),
            ),
            CoxeterDiagramType::E(rank @ 4..) => Self::new(
                rank,
                (2..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((0, 2), 3), ((1, 3), 3)]),
            ),
            CoxeterDiagramType::F(rank @ 3..) => Self::new(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((rank - 3, rank - 2), 4), ((rank - 2, rank - 1), 3)]),
            ),
            CoxeterDiagramType::G(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((rank - 2, rank - 1), 6)]),
            ),
            CoxeterDiagramType::H(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), 3))
                    .chain([((rank - 2, rank - 1), 5)]),
            ),
            CoxeterDiagramType::I2(label @ 3..) => Self::new(2, [((0, 1), label)]),
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

#[cfg(test)]
mod test {
    use super::{CoxeterDiagram, CoxeterDiagramType, CoxeterGroupInfo, CoxeterGroupType};

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

    #[test]
    fn test_coxeter_group_type() {
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::H(3)).coxeter_group_type(),
            CoxeterGroupType::Elliptic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::H(4)).coxeter_group_type(),
            CoxeterGroupType::Elliptic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::H(5)).coxeter_group_type(),
            CoxeterGroupType::Hyperbolic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::G(2)).coxeter_group_type(),
            CoxeterGroupType::Elliptic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::G(3)).coxeter_group_type(),
            CoxeterGroupType::Parabolic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::G(4)).coxeter_group_type(),
            CoxeterGroupType::Hyperbolic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(8)).coxeter_group_type(),
            CoxeterGroupType::Elliptic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(9)).coxeter_group_type(),
            CoxeterGroupType::Parabolic,
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(10)).coxeter_group_type(),
            CoxeterGroupType::Hyperbolic,
        );
    }
}

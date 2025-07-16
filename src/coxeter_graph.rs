use itertools::Itertools;

use super::alg::*;
use super::cyclotomic::Cyclotomic;
use super::square_matrix::SquareMatrix;

pub enum CoxeterDiagramNodeType {
    Inactive,
    Active,
}

#[derive(Clone)]
pub enum CoxeterGraphEdgeType {
    SelfLoop,
    Unconnected,
    Label(LabelType),
}

#[derive(Clone)]
pub enum LabelType {
    Implicit,
    Explicit(u32),
    ExplicitInfinite,
}

impl CoxeterGraphEdgeType {
    fn coxeter_matrix_entry(&self) -> Cyclotomic<i32, u32> {
        match self {
            &Self::SelfLoop => Cyclotomic::one().add(Cyclotomic::one()),
            &Self::Unconnected => Cyclotomic::zero(),
            &Self::Label(LabelType::Implicit) => Cyclotomic::one().neg(),
            &Self::Label(LabelType::Explicit(label)) => Cyclotomic::root_of_unity(2 * label, 1)
                .add(Cyclotomic::root_of_unity(2 * label, 2 * label - 1))
                .neg(),
            &Self::Label(LabelType::ExplicitInfinite) => {
                Cyclotomic::one().add(Cyclotomic::one()).neg()
            }
        }
    }
}

pub struct CoxeterGraph(petgraph::graph::UnGraph<(), LabelType>);

impl CoxeterGraph {
    pub fn from_edges<I>(rank: usize, edges: I) -> Self
    where
        I: IntoIterator<Item = ((usize, usize), LabelType)>,
    {
        let mut graph = petgraph::graph::UnGraph::new_undirected();
        let nodes = (0..rank).map(|_| graph.add_node(())).collect::<Vec<_>>();
        edges.into_iter().for_each(|((i, j), edge)| {
            graph.add_edge(nodes[i], nodes[j], edge);
        });
        Self(graph)
    }

    pub fn rank(&self) -> usize {
        self.0.node_count()
    }

    pub fn coxeter_matrix(&self) -> SquareMatrix<Cyclotomic<i32, u32>> {
        let indexing = |index| petgraph::visit::NodeIndexable::from_index(&self.0, index);
        SquareMatrix::from_fn(self.rank(), |i, j| {
            (if i == j {
                CoxeterGraphEdgeType::SelfLoop
            } else {
                self.0
                    .find_edge(indexing(i), indexing(j))
                    .map(|edge_index| self.0.edge_weight(edge_index))
                    .flatten()
                    .cloned()
                    .map(CoxeterGraphEdgeType::Label)
                    .unwrap_or(CoxeterGraphEdgeType::Unconnected)
            })
            .coxeter_matrix_entry()
        })
    }

    pub fn coxeter_element_embedding(&self) -> SquareMatrix<Cyclotomic<i32, u32>> {
        let coxeter_matrix = self.coxeter_matrix();
        let dimension = coxeter_matrix.dimension();
        let neg_strict_upper_triangular_matrix = SquareMatrix::from_fn(dimension, |i, j| {
            if i < j {
                coxeter_matrix.get(i, j).clone().neg()
            } else {
                Cyclotomic::zero()
            }
        });
        -(SquareMatrix::one(dimension) - neg_strict_upper_triangular_matrix.clone().transpose())
            * std::iter::successors(Some(SquareMatrix::one(dimension)), |matrix| {
                Some(neg_strict_upper_triangular_matrix.clone() * matrix.clone())
            })
            .take(dimension)
            .fold(SquareMatrix::zero(dimension), std::ops::Add::add)
    }

    pub fn coxeter_number(&self) -> u32 {
        println!(">>> {:#?}", self.coxeter_element_embedding());
        self.coxeter_element_embedding().order() as u32
    }

    pub fn exponents(&self) -> Vec<u32> {
        let coxeter_element_embedding = self.coxeter_element_embedding();
        let dimension = coxeter_element_embedding.dimension();
        let coxeter_number = coxeter_element_embedding.order() as u32;
        (0..coxeter_number)
            .filter(|&exponent| {
                println!("{exponent}");
                Cyclotomic::is_zero(
                    &(SquareMatrix::embed(
                        dimension,
                        Cyclotomic::root_of_unity(coxeter_number, exponent),
                    ) - coxeter_element_embedding.clone())
                    .determinant(),
                )
            })
            .collect()
    }

    pub fn degrees(&self) -> Vec<u32> {
        self.exponents()
            .into_iter()
            .map(|exponent| exponent + 1)
            .collect()
    }
}

pub enum CoxeterGroupType {
    A(usize),
    B(usize),
    C(usize),
    D(usize),
    E(usize),
    F(usize),
    G(usize),
    H(usize),
    I2(LabelType),
}

impl From<CoxeterGroupType> for CoxeterGraph {
    fn from(value: CoxeterGroupType) -> Self {
        match value {
            CoxeterGroupType::A(rank @ 1..) => CoxeterGraph::from_edges(
                rank,
                (0..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit)),
            ),
            CoxeterGroupType::B(rank @ 2..) => CoxeterGraph::from_edges(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([((rank - 2, rank - 1), LabelType::Explicit(4))]),
            ),
            CoxeterGroupType::C(rank @ 2..) => CoxeterGraph::from_edges(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([((rank - 2, rank - 1), LabelType::Explicit(4))]),
            ),
            CoxeterGroupType::D(rank @ 3..) => CoxeterGraph::from_edges(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([
                        ((rank - 3, rank - 2), LabelType::Implicit),
                        ((rank - 3, rank - 1), LabelType::Implicit),
                    ]),
            ),
            CoxeterGroupType::E(rank @ 4..=8) => CoxeterGraph::from_edges(
                rank,
                (2..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([((0, 2), LabelType::Implicit), ((1, 3), LabelType::Implicit)]),
            ),
            CoxeterGroupType::F(rank @ 3..=4) => CoxeterGraph::from_edges(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([
                        ((rank - 3, rank - 2), LabelType::Explicit(4)),
                        ((rank - 2, rank - 1), LabelType::Implicit),
                    ]),
            ),
            CoxeterGroupType::G(rank @ 2..=2) => CoxeterGraph::from_edges(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([((rank - 2, rank - 1), LabelType::Explicit(6))]),
            ),
            CoxeterGroupType::H(rank @ 2..=4) => CoxeterGraph::from_edges(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), LabelType::Implicit))
                    .chain([((rank - 2, rank - 1), LabelType::Explicit(5))]),
            ),
            CoxeterGroupType::I2(label) => CoxeterGraph::from_edges(2, [((0, 1), label)]),
            _ => panic!(),
        }
    }
}

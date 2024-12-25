use itertools::Itertools;

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
    Unlabelled,
    Label4,
    Label5,
    Label6,
    FiniteLabel(u64),
    InfiniteLabel,
}

impl CoxeterGraphEdgeType {
    fn coxeter_matrix_entry(&self) -> Cyclotomic<i64> {
        match self {
            &Self::SelfLoop => Some(1),
            &Self::Unconnected => Some(2),
            &Self::Unlabelled => Some(3),
            &Self::Label4 => Some(4),
            &Self::Label5 => Some(5),
            &Self::Label6 => Some(6),
            &Self::FiniteLabel(label) => Some(label),
            &Self::InfiniteLabel => None,
        }
        .map(|label| {
            -(Cyclotomic::<i64>::root_of_unity(2 * label, 1)
                + Cyclotomic::<i64>::root_of_unity(2 * label, 2 * label - 1))
        })
        .unwrap_or(num::Zero::zero())
    }
}

pub struct CoxeterGraph(petgraph::graph::UnGraph<(), CoxeterGraphEdgeType>);

impl CoxeterGraph {
    pub fn from_edges<I>(rank: usize, edges: I) -> Self
    where
        I: IntoIterator<Item = ((usize, usize), CoxeterGraphEdgeType)>,
    {
        let mut graph = petgraph::graph::UnGraph::new_undirected();
        let nodes = (0..rank).map(|_| graph.add_node(())).collect::<Vec<_>>();
        edges.into_iter().for_each(|((i, j), edge)| {
            graph.add_edge(nodes[i], nodes[j], edge);
        });
        Self(graph)
    }

    pub fn from_type(group_type: CoxeterGroupType) -> Self {
        coxeter_graph_from_type_extended(group_type, 0)
    }

    pub fn from_type_extended(group_type: CoxeterGroupType, extension: usize) -> Self {
        coxeter_graph_from_type_extended(group_type, extension)
    }

    pub fn rank(&self) -> usize {
        self.0.node_count()
    }

    pub fn coxeter_matrix(&self) -> SquareMatrix<Cyclotomic<i64>> {
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
                    .unwrap_or(CoxeterGraphEdgeType::Unconnected)
            })
            .coxeter_matrix_entry()
        })
    }

    pub fn coxeter_element_embedding(&self) -> SquareMatrix<Cyclotomic<i64>> {
        let mut coxeter_matrix = self.coxeter_matrix();
        let dimension = coxeter_matrix.dimension();
        let neg_strict_upper_triangular_matrix = SquareMatrix::from_fn(dimension, |i, j| {
            if i < j {
                -coxeter_matrix.take(i, j)
            } else {
                num::Zero::zero()
            }
        });
        -(SquareMatrix::one(dimension) - neg_strict_upper_triangular_matrix.clone().transpose())
            * neg_strict_upper_triangular_matrix.inverse_one_minus_nilpotent()
    }

    pub fn coxeter_number(&self) -> u64 {
        self.coxeter_element_embedding().order_unipotent() as u64
    }

    pub fn exponents(&self) -> Vec<u64> {
        let coxeter_element_embedding = self.coxeter_element_embedding();
        let coxeter_number = coxeter_element_embedding.order_unipotent() as u64;
        let characteristic_polynomial = coxeter_element_embedding.characteristic_polynomial();
        (0..coxeter_number)
            .filter(|&exponent| {
                num::Zero::is_zero(
                    &characteristic_polynomial
                        .clone()
                        .eval(Cyclotomic::<i64>::root_of_unity(coxeter_number, exponent)),
                )
            })
            .collect()
    }

    pub fn degrees(&self) -> Vec<u64> {
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
    E6,
    E7,
    E8,
    F4,
    G2,
    H3,
    H4,
    I2(u64),
}

impl CoxeterGroupType {
    pub fn rank(&self) -> usize {
        match *self {
            Self::A(rank) => rank,
            Self::B(rank) => rank,
            Self::C(rank) => rank,
            Self::D(rank) => rank,
            Self::E6 => 6,
            Self::E7 => 7,
            Self::E8 => 8,
            Self::F4 => 4,
            Self::G2 => 2,
            Self::H3 => 3,
            Self::H4 => 4,
            Self::I2(_) => 2,
        }
    }
}

// https://en.wikipedia.org/wiki/Coxeter%E2%80%93Dynkin_diagram
fn coxeter_graph_from_type_extended(
    group_type: CoxeterGroupType,
    extension: usize,
) -> CoxeterGraph {
    let rank = group_type.rank();
    let (edges, joints) = match group_type {
        CoxeterGroupType::A(1) => (
            [].to_vec(),
            [((0, rank), CoxeterGraphEdgeType::InfiniteLabel)].to_vec(),
        ),
        CoxeterGroupType::A(rank) => (
            (0..rank)
                .tuple_windows()
                .map(|(i, j)| ((i, j), CoxeterGraphEdgeType::Unlabelled))
                .collect(),
            [
                ((0, rank), CoxeterGraphEdgeType::Unlabelled),
                ((rank - 1, rank), CoxeterGraphEdgeType::Unlabelled),
            ]
            .to_vec(),
        ),
        CoxeterGroupType::B(rank) if rank >= 3 => (
            (0..rank - 2)
                .tuple_windows()
                .map(|(i, j)| ((i, j), CoxeterGraphEdgeType::Unlabelled))
                .chain([
                    ((rank - 3, rank - 2), CoxeterGraphEdgeType::Label4),
                    ((0, rank - 1), CoxeterGraphEdgeType::Unlabelled),
                ])
                .collect(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::C(rank) if rank >= 2 => (
            (0..rank - 1)
                .tuple_windows()
                .map(|(i, j)| ((i, j), CoxeterGraphEdgeType::Unlabelled))
                .chain([((rank - 2, rank - 1), CoxeterGraphEdgeType::Label4)])
                .collect(),
            [((0, rank), CoxeterGraphEdgeType::Label4)].to_vec(),
        ),
        CoxeterGroupType::D(rank) if rank >= 4 => (
            (0..rank - 2)
                .tuple_windows()
                .map(|(i, j)| ((i, j), CoxeterGraphEdgeType::Unlabelled))
                .chain([
                    ((rank - 4, rank - 2), CoxeterGraphEdgeType::Unlabelled),
                    ((0, rank - 1), CoxeterGraphEdgeType::Unlabelled),
                ])
                .collect(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::E6 => (
            [
                ((0, 1), CoxeterGraphEdgeType::Unlabelled),
                ((1, 2), CoxeterGraphEdgeType::Unlabelled),
                ((2, 3), CoxeterGraphEdgeType::Unlabelled),
                ((1, 4), CoxeterGraphEdgeType::Unlabelled),
                ((4, 5), CoxeterGraphEdgeType::Unlabelled),
            ]
            .to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::E7 => (
            [
                ((0, 1), CoxeterGraphEdgeType::Unlabelled),
                ((1, 2), CoxeterGraphEdgeType::Unlabelled),
                ((2, 3), CoxeterGraphEdgeType::Unlabelled),
                ((3, 4), CoxeterGraphEdgeType::Unlabelled),
                ((4, 5), CoxeterGraphEdgeType::Unlabelled),
                ((2, 6), CoxeterGraphEdgeType::Unlabelled),
            ]
            .to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::E8 => (
            [
                ((0, 1), CoxeterGraphEdgeType::Unlabelled),
                ((1, 2), CoxeterGraphEdgeType::Unlabelled),
                ((2, 3), CoxeterGraphEdgeType::Unlabelled),
                ((3, 4), CoxeterGraphEdgeType::Unlabelled),
                ((4, 5), CoxeterGraphEdgeType::Unlabelled),
                ((5, 6), CoxeterGraphEdgeType::Unlabelled),
                ((4, 7), CoxeterGraphEdgeType::Unlabelled),
            ]
            .to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::F4 => (
            [
                ((0, 1), CoxeterGraphEdgeType::Unlabelled),
                ((1, 2), CoxeterGraphEdgeType::Label4),
                ((2, 3), CoxeterGraphEdgeType::Unlabelled),
            ]
            .to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::G2 => (
            [((0, 1), CoxeterGraphEdgeType::Label6)].to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::H3 => (
            [
                ((0, 1), CoxeterGraphEdgeType::Unlabelled),
                ((1, 2), CoxeterGraphEdgeType::Label5),
            ]
            .to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::H4 => (
            [
                ((0, 1), CoxeterGraphEdgeType::Unlabelled),
                ((1, 2), CoxeterGraphEdgeType::Unlabelled),
                ((2, 3), CoxeterGraphEdgeType::Label5),
            ]
            .to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(),
        ),
        CoxeterGroupType::I2(label) if label >= 3 => (
            [((0, 1), CoxeterGraphEdgeType::FiniteLabel(label))].to_vec(),
            [((0, rank), CoxeterGraphEdgeType::Unlabelled)].to_vec(), // check
        ),
        _ => panic!(),
    };
    if extension == 0 {
        CoxeterGraph::from_edges(rank, edges.into_iter())
    } else {
        CoxeterGraph::from_edges(
            rank + extension,
            edges.into_iter().chain(joints.into_iter()).chain(
                (rank + 1..extension)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), CoxeterGraphEdgeType::Unlabelled)),
            ),
        )
    }
}

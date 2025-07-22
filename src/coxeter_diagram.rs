use itertools::Itertools;

// use super::alg::*;
use super::cyclotomic::Cyclotomic;
use super::dynkin_diagram::DynkinDiagramType;
use super::square_matrix::SquareMatrix;

pub enum Label {
    Implicit,
    Explicit(u32),
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

    fn get_edge(&self, (i, j): (usize, usize)) -> Option<&Label> {
        self.0
            .find_edge(
                petgraph::visit::NodeIndexable::from_index(&self.0, i),
                petgraph::visit::NodeIndexable::from_index(&self.0, j),
            )
            .map(|edge_index| self.0.edge_weight(edge_index).unwrap())
    }

    pub fn coxeter_matrix(&self) -> SquareMatrix<u32> {
        SquareMatrix::from_fn(self.rank(), |(i, j)| {
            self.get_edge((i, j))
                .map(|label| match label {
                    Label::Implicit => 3,
                    Label::Explicit(label) => *label,
                })
                .unwrap_or_else(|| if i == j { 1 } else { 2 })
        })
    }

    pub fn schlafli_matrix(&self) -> SquareMatrix<Cyclotomic> {
        SquareMatrix::from_fn(self.rank(), |(i, j)| {
            self.get_edge((i, j))
                .map(|label| match label {
                    Label::Implicit => Cyclotomic::embed_scalar(-1),
                    Label::Explicit(label) => {
                        -(Cyclotomic::root_of_unity(2 * label, 1)
                            + Cyclotomic::root_of_unity(2 * label, 2 * label - 1))
                    }
                })
                .unwrap_or_else(|| {
                    if i == j {
                        Cyclotomic::embed_scalar(2)
                    } else {
                        Cyclotomic::embed_scalar(0)
                    }
                })
        })
    }

    pub fn coxeter_element(&self) -> SquareMatrix<Cyclotomic> {
        let schlafli_matrix = self.schlafli_matrix();
        let dimension = schlafli_matrix.dimension();
        let neg_strict_upper_triangular_matrix = SquareMatrix::from_fn(dimension, |(i, j)| {
            if i < j {
                -schlafli_matrix.get((i, j)).clone()
            } else {
                Cyclotomic::embed_scalar(0)
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
        self.coxeter_element().order() as u32
    }

    pub fn exponents(&self) -> Vec<u32> {
        let coxeter_element = self.coxeter_element();
        let dimension = coxeter_element.dimension();
        let coxeter_number = coxeter_element.order() as u32;
        (0..coxeter_number)
            .filter(|&exponent| {
                Cyclotomic::is_zero(
                    &(SquareMatrix::embed(
                        dimension,
                        Cyclotomic::root_of_unity(coxeter_number, exponent),
                    ) - coxeter_element.clone())
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

pub enum CoxeterDiagramType {
    A(usize),
    B(usize),
    C(usize),
    D(usize),
    E(usize),
    F(usize),
    G(usize),
    H(usize),
    I2(u32),
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

use super::square_matrix::SquareMatrixRingBase;
use feanor_math::{matrix::OwnedMatrix, ring::RingBase};
use itertools::Itertools;

pub enum Bond {
    Single,
    Double,
    Triple,
}

pub enum DynkinDiagramType {
    A(usize),
    B(usize),
    C(usize),
    D(usize),
    E(usize),
    F(usize),
    G(usize),
}

pub struct DynkinDiagram(petgraph::graph::DiGraph<(), Bond>);

impl DynkinDiagram {
    pub fn new<I>(rank: usize, bonds: I) -> Self
    where
        I: IntoIterator<Item = ((usize, usize), Bond)>,
    {
        let mut graph = petgraph::graph::DiGraph::new();
        (0..rank).for_each(|_| {
            graph.add_node(());
        });
        bonds.into_iter().for_each(|((i, j), bond)| {
            graph.update_edge(
                petgraph::graph::NodeIndex::new(i),
                petgraph::graph::NodeIndex::new(j),
                bond,
            );
        });
        Self(graph)
    }

    fn get_edge(&self, i: usize, j: usize) -> Option<&Bond> {
        self.0
            .find_edge(
                petgraph::visit::NodeIndexable::from_index(&self.0, i),
                petgraph::visit::NodeIndexable::from_index(&self.0, j),
            )
            .map(|edge_index| self.0.edge_weight(edge_index).unwrap())
    }

    pub fn rank(&self) -> usize {
        self.0.node_count()
    }

    pub fn cartan_matrix(
        &self,
    ) -> <SquareMatrixRingBase<feanor_math::primitive_int::StaticRing<i64>> as RingBase>::Element
    {
        let rank = self.rank();
        OwnedMatrix::from_fn(rank, rank, |i, j| {
            self.get_edge(i, j)
                .map(|bond| match bond {
                    Bond::Single => -1,
                    Bond::Double => -1,
                    Bond::Triple => -1,
                })
                .or_else(|| {
                    self.get_edge(j, i).map(|bond| match bond {
                        Bond::Single => -1,
                        Bond::Double => -2,
                        Bond::Triple => -3,
                    })
                })
                .unwrap_or_else(|| if i == j { 2 } else { 0 })
        })
    }
}

impl From<DynkinDiagramType> for DynkinDiagram {
    fn from(value: DynkinDiagramType) -> Self {
        match value {
            DynkinDiagramType::A(rank @ 1..) => Self::new(
                rank,
                (0..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single)),
            ),
            DynkinDiagramType::B(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single))
                    .chain([((rank - 1, rank - 2), Bond::Double)]),
            ),
            DynkinDiagramType::C(rank @ 2..) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single))
                    .chain([((rank - 2, rank - 1), Bond::Double)]),
            ),
            DynkinDiagramType::D(rank @ 3..) => Self::new(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single))
                    .chain([
                        ((rank - 3, rank - 2), Bond::Single),
                        ((rank - 3, rank - 1), Bond::Single),
                    ]),
            ),
            DynkinDiagramType::E(rank @ 4..=8) => Self::new(
                rank,
                (2..rank)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single))
                    .chain([((0, 2), Bond::Single), ((1, 3), Bond::Single)]),
            ),
            DynkinDiagramType::F(rank @ 3..=4) => Self::new(
                rank,
                (0..rank - 2)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single))
                    .chain([
                        ((rank - 2, rank - 3), Bond::Double),
                        ((rank - 2, rank - 1), Bond::Single),
                    ]),
            ),
            DynkinDiagramType::G(rank @ 2..=2) => Self::new(
                rank,
                (0..rank - 1)
                    .tuple_windows()
                    .map(|(i, j)| ((i, j), Bond::Single))
                    .chain([((rank - 2, rank - 1), Bond::Triple)]),
            ),
            _ => panic!(),
        }
    }
}

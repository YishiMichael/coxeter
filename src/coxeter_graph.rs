pub enum CoxeterGraphNodeType {
    Inactive,
    Active,
}

pub struct CoxeterGraph<T>(petgraph::graph::UnGraph<(), T>);

impl<T> CoxeterGraph<T> {
    pub fn from_edges<I>(iterable: I) -> Self
    where
        I: IntoIterator<Item = ((u32, u32), T)>,
    {
        Self(petgraph::graph::UnGraph::<(), T>::from_edges(
            iterable.into_iter().map(|((a, b), weight)| (a, b, weight)),
        ))
    }
}

pub mod alg;
pub mod big_rational;
pub mod big_uint;
pub mod coxeter_graph;
pub mod cyclotomic;
pub mod graded_algebra;
pub mod polynomial;
pub mod prufer_group;
pub mod square_matrix;

#[cfg(test)]
mod tests {
    use super::alg;
    use super::alg::AdditiveIdentity;
    use super::alg::AdditiveInverse;

    #[test]
    fn test_cyclotomic() {
        use super::cyclotomic::Cyclotomic;
        assert!(Cyclotomic::root_of_unity(3, 0) == Cyclotomic::root_of_unity(2, 0));
        assert!(alg::AdditiveMagma::sum(
            [
                Cyclotomic::root_of_unity(3, 1),
                Cyclotomic::root_of_unity(3, 2),
                Cyclotomic::root_of_unity(3, 0),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(alg::AdditiveMagma::sum(
            [
                Cyclotomic::root_of_unity(2, 0),
                Cyclotomic::root_of_unity(6, 2),
                Cyclotomic::root_of_unity(3, 2),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(alg::AdditiveMagma::sum(
            [
                Cyclotomic::root_of_unity(5, 1),
                Cyclotomic::root_of_unity(5, 2),
                Cyclotomic::root_of_unity(5, 3),
                Cyclotomic::root_of_unity(5, 4),
                Cyclotomic::root_of_unity(3, 1).neg(),
                Cyclotomic::root_of_unity(3, 2).neg(),
            ]
            .into_iter()
        )
        .is_zero());
    }

    #[test]
    fn test_coxeter_graph() {
        use super::coxeter_graph::CoxeterGraph;
        use super::coxeter_graph::CoxeterGroupType;
        assert_eq!(
            CoxeterGraph::from_type(CoxeterGroupType::B(3)).coxeter_number(),
            6
        );
        assert_eq!(
            CoxeterGraph::from_type(CoxeterGroupType::H3).coxeter_number(),
            10
        );
        assert_eq!(
            CoxeterGraph::from_type(CoxeterGroupType::E6).coxeter_number(),
            12
        );
        assert_eq!(
            CoxeterGraph::from_type(CoxeterGroupType::E7).coxeter_number(),
            18
        );
        assert_eq!(
            CoxeterGraph::from_type(CoxeterGroupType::E8).coxeter_number(),
            30
        );
    }
}

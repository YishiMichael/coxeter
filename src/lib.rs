pub mod alg;
pub mod coxeter_graph;
pub mod cyclotomic;
pub mod graded_algebra;
pub mod prufer_group;
pub mod square_matrix;

#[cfg(test)]
mod tests {
    use super::alg::*;

    #[test]
    fn test_cyclotomic() {
        use super::cyclotomic::Cyclotomic;
        assert!(Cyclotomic::<i32, u32>::root_of_unity(3, 0) == Cyclotomic::root_of_unity(2, 0));
        assert!(Cyclotomic::sum(
            [
                Cyclotomic::<i32, u32>::root_of_unity(3, 1),
                Cyclotomic::<i32, u32>::root_of_unity(3, 2),
                Cyclotomic::<i32, u32>::root_of_unity(3, 0),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(Cyclotomic::sum(
            [
                Cyclotomic::<i32, u32>::root_of_unity(2, 0),
                Cyclotomic::<i32, u32>::root_of_unity(6, 2),
                Cyclotomic::<i32, u32>::root_of_unity(3, 2),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(Cyclotomic::sum(
            [
                Cyclotomic::<i32, u32>::root_of_unity(5, 1),
                Cyclotomic::<i32, u32>::root_of_unity(5, 2),
                Cyclotomic::<i32, u32>::root_of_unity(5, 3),
                Cyclotomic::<i32, u32>::root_of_unity(5, 4),
                Cyclotomic::<i32, u32>::root_of_unity(3, 1).neg(),
                Cyclotomic::<i32, u32>::root_of_unity(3, 2).neg(),
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
            CoxeterGraph::from(CoxeterGroupType::B(3)).coxeter_number(),
            6
        );
        assert_eq!(
            CoxeterGraph::from(CoxeterGroupType::H(3)).coxeter_number(),
            10
        );
        assert_eq!(
            CoxeterGraph::from(CoxeterGroupType::E(6)).coxeter_number(),
            12
        );
        assert_eq!(
            CoxeterGraph::from(CoxeterGroupType::E(7)).coxeter_number(),
            18
        );
        assert_eq!(
            CoxeterGraph::from(CoxeterGroupType::E(8)).coxeter_number(),
            30
        );
    }

    #[test]
    fn test_coxeter_graph_degrees() {
        use super::coxeter_graph::CoxeterGraph;
        use super::coxeter_graph::CoxeterGroupType;
        assert_eq!(
            CoxeterGraph::from(CoxeterGroupType::E(8)).degrees(),
            &[2, 8, 12, 14, 18, 20, 24, 30]
        );
    }
}

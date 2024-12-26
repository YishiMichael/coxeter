pub mod alg;
pub mod big_rational;
pub mod big_uint;
// pub mod coxeter_graph;
pub mod cyclotomic;
pub mod graded_algebra;
pub mod polynomial;
pub mod prufer_group;
pub mod square_matrix;

#[cfg(test)]
mod tests {

    #[test]
    fn test_cyclotomic() {
        // use super::graded_algebra::GradedAlgebraMeta;

        // use super::cyclotomic::Cyclotomic;
        // use super::cyclotomic::CyclotomicMeta;
        // assert!(Cyclotomic::<i64>::root_of_unity(3, 0) == Cyclotomic::<i64>::root_of_unity(2, 0));
        // dbg!(CyclotomicMeta::annihilate(
        //     Cyclotomic::<i64>::root_of_unity(3, 0)
        //         + Cyclotomic::<i64>::root_of_unity(3, 1)
        //         + Cyclotomic::<i64>::root_of_unity(3, 2)
        // ));
        // assert!(
        //     -Cyclotomic::<i64>::root_of_unity(3, 0)
        //         == Cyclotomic::<i64>::root_of_unity(3, 1) + Cyclotomic::<i64>::root_of_unity(3, 2)
        // );
        // assert!(
        //     -Cyclotomic::<i64>::root_of_unity(2, 0)
        //         == Cyclotomic::<i64>::root_of_unity(6, 2) + Cyclotomic::<i64>::root_of_unity(3, 2)
        // );
        // dbg!(CyclotomicMeta::annihilate(
        //     Cyclotomic::<i64>::root_of_unity(5, 1)
        //         + Cyclotomic::<i64>::root_of_unity(5, 2)
        //         + Cyclotomic::<i64>::root_of_unity(5, 3)
        //         + Cyclotomic::<i64>::root_of_unity(5, 4)
        //         - Cyclotomic::<i64>::root_of_unity(3, 1)
        //         - Cyclotomic::<i64>::root_of_unity(3, 2)
        // ));
        // assert!(
        //     Cyclotomic::<i64>::root_of_unity(5, 1)
        //         + Cyclotomic::<i64>::root_of_unity(5, 2)
        //         + Cyclotomic::<i64>::root_of_unity(5, 3)
        //         + Cyclotomic::<i64>::root_of_unity(5, 4)
        //         == Cyclotomic::<i64>::root_of_unity(3, 1) + Cyclotomic::<i64>::root_of_unity(3, 2)
        // );
    }

    #[test]
    fn test_coxeter_graph() {
        // use super::coxeter_graph::CoxeterGraph;
        // use super::coxeter_graph::CoxeterGroupType;
        // assert_eq!(
        //     CoxeterGraph::from_type(CoxeterGroupType::B(3)).coxeter_number(),
        //     6
        // );
        // assert_eq!(
        //     CoxeterGraph::from_type(CoxeterGroupType::H3).coxeter_number(),
        //     10
        // );
        // assert_eq!(
        //     CoxeterGraph::from_type(CoxeterGroupType::E6).coxeter_number(),
        //     12
        // );
        // assert_eq!(
        //     CoxeterGraph::from_type(CoxeterGroupType::E7).coxeter_number(),
        //     18
        // );
        // assert_eq!(
        //     CoxeterGraph::from_type(CoxeterGroupType::E8).coxeter_number(),
        //     30
        // );
    }
}

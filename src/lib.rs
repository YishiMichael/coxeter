pub mod coxeter_diagram;
pub mod cyclotomic;
pub mod dynkin_diagram;
pub mod square_matrix;

#[cfg(test)]
mod tests {
    // use super::alg::*;

    #[test]
    fn test_cyclotomic() {
        use super::cyclotomic::Cyclotomic;
        assert!(Cyclotomic::root_of_unity(3, 0) == Cyclotomic::root_of_unity(2, 0));
        assert!(Cyclotomic::sum(
            [
                Cyclotomic::root_of_unity(3, 1),
                Cyclotomic::root_of_unity(3, 2),
                Cyclotomic::root_of_unity(3, 0),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(Cyclotomic::sum(
            [
                Cyclotomic::root_of_unity(3, 1),
                Cyclotomic::root_of_unity(9, 6),
                Cyclotomic::root_of_unity(3, 0),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(Cyclotomic::sum(
            [
                Cyclotomic::root_of_unity(2, 0),
                Cyclotomic::root_of_unity(6, 2),
                Cyclotomic::root_of_unity(3, 2),
            ]
            .into_iter()
        )
        .is_zero());
        assert!(Cyclotomic::sum(
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
    fn test_coxeter_numbers() {
        use super::coxeter_diagram::CoxeterDiagram;
        use super::coxeter_diagram::CoxeterDiagramType;
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::B(3)).coxeter_number(),
            6
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::H(3)).coxeter_number(),
            10
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::G(2)).coxeter_number(),
            6
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::F(4)).coxeter_number(),
            12
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(6)).coxeter_number(),
            12
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(7)).coxeter_number(),
            18
        );
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(8)).coxeter_number(),
            30
        );
    }

    #[test]
    fn test_g2_degrees() {
        use super::coxeter_diagram::CoxeterDiagram;
        use super::coxeter_diagram::CoxeterDiagramType;
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::G(2)).degrees(),
            &[2, 6]
        );
    }

    #[test]
    fn test_f4_degrees() {
        use super::coxeter_diagram::CoxeterDiagram;
        use super::coxeter_diagram::CoxeterDiagramType;
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::F(4)).degrees(),
            &[2, 6, 8, 12]
        );
    }

    #[test]
    fn test_e6_degrees() {
        use super::coxeter_diagram::CoxeterDiagram;
        use super::coxeter_diagram::CoxeterDiagramType;
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(6)).degrees(),
            &[2, 5, 6, 8, 9, 12]
        );
    }

    #[test]
    fn test_e7_degrees() {
        use super::coxeter_diagram::CoxeterDiagram;
        use super::coxeter_diagram::CoxeterDiagramType;
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(7)).degrees(),
            &[2, 6, 8, 10, 12, 14, 18]
        );
    }

    #[test]
    fn test_e8_degrees() {
        use super::coxeter_diagram::CoxeterDiagram;
        use super::coxeter_diagram::CoxeterDiagramType;
        assert_eq!(
            CoxeterDiagram::from(CoxeterDiagramType::E(8)).degrees(),
            &[2, 8, 12, 14, 18, 20, 24, 30]
        );
    }
}

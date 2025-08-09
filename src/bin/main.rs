use coxeter::*;

fn main() {
    for (rank, depth_bound) in [
        (1, 6),
        (2, 6),
        (3, 3),
        (4, 2),
        (5, 2),
        (6, 2),
        (7, 2),
        (8, 2),
        (9, 2),
        (10, 2),
    ] {
        for (coxeter_diagram, coxeter_group_type) in
            coxeter_diagram::CoxeterDiagram::enumerate_trees_classified(rank, depth_bound)
        {
            if coxeter_group_type != coxeter_diagram::CoxeterGroupType::Hyperbolic {
                dbg!(coxeter_diagram.rank());
                dbg!(coxeter_diagram);
                dbg!(coxeter_group_type);
            }
        }
    }
}

use coxeter::*;

fn main() {
    let mut verified_non_elliptic_coxeter_diagrams = Vec::<coxeter_diagram::CoxeterDiagram>::new();
    for (rank, rank_iter) in coxeter_diagram::CoxeterDiagram::enumerate_trees_by_rank_and_depth() {
        for (depth, coxeter_diagrams) in rank_iter {
            let mut has_elliptic = false;
            for coxeter_diagram in coxeter_diagrams {
                if verified_non_elliptic_coxeter_diagrams
                    .iter()
                    .any(|non_elliptic| non_elliptic.is_subgraph_of(&coxeter_diagram))
                {
                    continue;
                }
                let coxeter_group_type = coxeter_diagram.coxeter_group_type();
                if coxeter_group_type == coxeter_diagram::CoxeterGroupType::Elliptic
                    || coxeter_group_type == coxeter_diagram::CoxeterGroupType::Parabolic
                {
                    println!("{:?} -> {:?}", &coxeter_diagram, coxeter_group_type);
                }
                if coxeter_group_type == coxeter_diagram::CoxeterGroupType::Elliptic {
                    has_elliptic = true;
                } else {
                    verified_non_elliptic_coxeter_diagrams.push(coxeter_diagram);
                }
            }
            if !has_elliptic || depth > 12 {
                break;
            }
        }
        if rank > 12 {
            break;
        }
    }
}

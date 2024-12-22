use alga::general as alg;

// natural number
pub trait NN: num::Integer + num::Unsigned + Clone {}

fn simplify<N: NN>(numerator: N, denominator: N) -> (N, N) {
    let gcd = numerator.gcd(&denominator);
    (numerator / gcd.clone(), denominator / gcd)
}

// Convention: 0 is represented as 0 / 1
#[derive(Clone, PartialEq)]
pub struct PruferGroup<N: NN> {
    numerator: N,
    denominator: N,
}

impl<N: NN> PruferGroup<N> {
    pub fn new(numerator: N, denominator: N) -> Self {
        Self {
            numerator,
            denominator,
        }
    }
}

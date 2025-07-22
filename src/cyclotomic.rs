use itertools::Itertools;

#[derive(Clone)]
pub struct Cyclotomic(ndarray::Array1<i32>);

impl Cyclotomic {
    pub fn root_of_unity(order: u32, exponent: u32) -> Self {
        assert!(order != 0);
        let mut coefficients = ndarray::Array1::zeros((order as usize,));
        coefficients[exponent as usize % order as usize] = 1;
        Self(coefficients)
        // let factorization: Vec<_> = PrimesIter::default()
        //     .scan(order, |n, p| {
        //         (*n != 1).then(|| {
        //             (
        //                 p,
        //                 std::iter::from_fn(|| (*n % p == 0).then_some(()).inspect(|_| *n /= p))
        //                     .count() as u32,
        //             )
        //         })
        //     })
        //     .map(|(p, d)| {
        //         let q = p.pow(d);
        //         let q_numer = Self::modpow(order / q, phi_q - 1, q) * exponent % q;

        //         // d.checked_sub(1).map(|d_m1| {
        //         //     let q = p.pow(d);
        //         //     let phi_q = (p - 1) * p.pow(d_m1);
        //         //     let p_digits: Vec<_> = (0..d)
        //         //         .scan(Self::modpow(order / q, phi_q - 1, q) * exponent, |n, _| {
        //         //             let rem = *n % p;
        //         //             *n /= p;
        //         //             Some(rem)
        //         //         })
        //         //         .collect();
        //         //     let slice_info_elem = if *p_digits.first().unwrap() == 0 {
        //         //         ndarray::SliceInfoElem::Slice { start: (), end: (), step: 1 }
        //         //     } else {
        //         //         ndarray::SliceInfoElem::Index(
        //         //             (p_digits.first().unwrap() - 1
        //         //                 + p_digits
        //         //                     .iter()
        //         //                     .skip(1)
        //         //                     .rev()
        //         //                     .fold(0, |acc, p_digit| acc * p + *p_digit)
        //         //                     * (p - 1)) as isize,
        //         //         )
        //         //     };
        //         // })
        //     })
        //     .collect();
        // // let shape =
    }

    fn modpow(base: u32, exp: u32, modulus: u32) -> u32 {
        std::iter::successors(Some(exp), |&e| (e != 0).then(|| e / 2))
            .fold((base, 1), |(base, result), e| {
                (
                    base * base % modulus,
                    if e % 2 == 1 {
                        result * base % modulus
                    } else {
                        result
                    },
                )
            })
            .1
    }

    fn reverse(n: u32, base: u32, exp: u32) -> u32 {
        (0..exp)
            .scan(n, |n, _| {
                let rem = *n % base;
                *n /= base;
                Some(rem)
            })
            .collect_vec()
            .into_iter()
            .rev()
            .fold(0, |n, rem| n * base + rem)
    }
}

#[derive(Default)]
struct PrimesIter(Vec<u32>);

impl Iterator for PrimesIter {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        let prime = (self.0.last().cloned().unwrap_or(1) + 1..)
            .find_or_last(|&n| {
                self.0
                    .iter()
                    .take_while(|&&p| p * p <= n)
                    .all(|&p| n % p != 0)
            })
            .unwrap();
        self.0.push(prime);
        Some(prime)
    }
}

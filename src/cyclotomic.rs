use itertools::Itertools;
use ndarray::{s, Array1};
use num::Integer;

#[derive(Clone)]
pub struct Cyclotomic(Array1<i32>);

impl Cyclotomic {
    pub fn embed_scalar(scalar: i32) -> Self {
        Self(scalar * Array1::ones((1,)))
    }

    pub fn root_of_unity(order: u32, exponent: u32) -> Self {
        assert!(order != 0);
        let mut coeffs = Array1::zeros((order as usize,));
        coeffs[exponent as usize % order as usize] = 1;
        Self(coeffs)
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

    pub fn is_zero(&self) -> bool {
        let order = self.0.len();
        (2_usize..)
            .scan(order, |n, p| {
                (*n != 1).then(|| {
                    std::iter::from_fn(|| (*n % p == 0).then(|| *n /= p))
                        .next()
                        .map(|()| p)
                })
            })
            .filter_map(|option_p| option_p)
            .all(|p| {
                self.0
                    .to_shape((order / p, p))
                    .unwrap()
                    .axis_iter(ndarray::Axis(0))
                    .all(|slice| slice.sum() == 0)
            })
    }

    // fn modpow(base: u32, exp: u32, modulus: u32) -> u32 {
    //     std::iter::successors(Some(exp), |&e| (e != 0).then(|| e / 2))
    //         .fold((base, 1), |(base, result), e| {
    //             (
    //                 base * base % modulus,
    //                 if e % 2 == 1 {
    //                     result * base % modulus
    //                 } else {
    //                     result
    //                 },
    //             )
    //         })
    //         .1
    // }

    // fn reverse(n: u32, base: u32, exp: u32) -> u32 {
    //     (0..exp)
    //         .scan(n, |n, _| {
    //             let rem = *n % base;
    //             *n /= base;
    //             Some(rem)
    //         })
    //         .collect_vec()
    //         .into_iter()
    //         .rev()
    //         .fold(0, |n, rem| n * base + rem)
    // }
}

impl PartialEq for Cyclotomic {
    fn eq(&self, other: &Self) -> bool {
        (self - other).is_zero()
    }
}

impl Eq for Cyclotomic {}

impl std::ops::Add for &Cyclotomic {
    type Output = Cyclotomic;

    fn add(self, rhs: &Cyclotomic) -> Self::Output {
        let order = self.0.len().lcm(&rhs.0.len());
        let mut coeffs_lhs = Array1::zeros((order,));
        coeffs_lhs
            .slice_mut(s![..;order / self.0.len()])
            .assign(&self.0);
        let mut coeffs_rhs = Array1::zeros((order,));
        coeffs_rhs
            .slice_mut(s![..;order / rhs.0.len()])
            .assign(&rhs.0);
        Cyclotomic(coeffs_lhs + coeffs_rhs)
    }
}

impl std::ops::Sub for &Cyclotomic {
    type Output = Cyclotomic;

    fn sub(self, rhs: &Cyclotomic) -> Self::Output {
        let order = self.0.len().lcm(&rhs.0.len());
        let mut coeffs_lhs = Array1::zeros((order,));
        coeffs_lhs
            .slice_mut(s![..;order / self.0.len()])
            .assign(&self.0);
        let mut coeffs_rhs = Array1::zeros((order,));
        coeffs_rhs
            .slice_mut(s![..;order / rhs.0.len()])
            .assign(&rhs.0);
        Cyclotomic(coeffs_lhs - coeffs_rhs)
    }
}

impl std::ops::Mul for &Cyclotomic {
    type Output = Cyclotomic;

    fn mul(self, rhs: &Cyclotomic) -> Self::Output {
        let order = self.0.len().lcm(&rhs.0.len());
        let step_lhs = order / self.0.len();
        let step_rhs = order / rhs.0.len();
        let mut coeffs = Array1::zeros((order,));
        self.0
            .iter()
            .enumerate()
            .filter(|(_, &coeff_lhs)| coeff_lhs != 0)
            .cartesian_product(
                rhs.0
                    .iter()
                    .enumerate()
                    .filter(|(_, &coeff_rhs)| coeff_rhs != 0),
            )
            .for_each(|((index_lhs, &coeff_lhs), (index_rhs, &coeff_rhs))| {
                coeffs[(index_lhs * step_lhs + index_rhs * step_rhs) % order] +=
                    coeff_lhs * coeff_rhs
            });
        Cyclotomic(coeffs)
    }
}

impl std::ops::Add for Cyclotomic {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl std::ops::Sub for Cyclotomic {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl std::ops::Mul for Cyclotomic {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl std::ops::Neg for Cyclotomic {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

// #[derive(Default)]
// struct PrimesIter(Vec<u32>);

// impl Iterator for PrimesIter {
//     type Item = u32;

//     fn next(&mut self) -> Option<Self::Item> {
//         let prime = (self.0.last().cloned().unwrap_or(1) + 1..)
//             .find_or_last(|&n| {
//                 self.0
//                     .iter()
//                     .take_while(|&&p| p * p <= n)
//                     .all(|&p| n % p != 0)
//             })
//             .unwrap();
//         self.0.push(prime);
//         Some(prime)
//     }
// }

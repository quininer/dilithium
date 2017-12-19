extern crate hex;

use self::hex::FromHexError;
use super::*;
use byteorder::{ ByteOrder, BigEndian };
use itertools::Itertools;
use polyvec::{ PolyVecL, PolyVecK };
use params::{
    N, K, L,
    SEEDBYTES, CRHBYTES
};

const TEST_VECTORS: &str = include_str!("../tests/testvectors.txt");

struct TestVector {
    seed: ([u8; SEEDBYTES], [u8; CRHBYTES]),
    mat: [PolyVecL; K],
    s: PolyVecL,
    y: PolyVecL,
    w1: PolyVecK,
    c: [u32; N]
}

impl Default for TestVector {
    fn default() -> Self {
        TestVector {
            seed: ([0; SEEDBYTES], [0; CRHBYTES]),
            mat: [PolyVecL::default(); K],
            s: PolyVecL::default(),
            y: PolyVecL::default(),
            w1: PolyVecK::default(),
            c: [0; N]
        }
    }
}

fn parse_testvectors() -> Result<Vec<TestVector>, FromHexError> {
    let mut testvectors = Vec::new();

    for testvector in TEST_VECTORS.lines().chunks(7).into_iter() {
        let mut tv = TestVector::default();

        for (key, val) in testvector
            .map(|line| line.split('='))
            .map(|mut split| (split.next(), split.last()))
            .filter_map(|(key, val)| key.and_then(|key| val.map(|val| (key.trim(), val.trim()))))
        {
            match key {
                "count" => (),
                "seed" => {
                    let seed = hex::decode(val)?;
                    let (rho, mu) = seed.split_at(SEEDBYTES);
                    tv.seed.0.copy_from_slice(rho);
                    tv.seed.1.copy_from_slice(mu);
                },
                "mat" => {
                    let mut mat = [0; K * L * N];
                    let mut i = 0;
                    BigEndian::read_u32_into(&hex::decode(val)?, &mut mat);
                    for j in 0..K {
                        for k in 0..L {
                            for l in 0..N {
                                tv.mat[j][k][l] = mat[i];
                                i += 1;
                            }
                        }
                    }
                },
                "s" => {
                    let mut s = [0; L * N];
                    let mut i = 0;
                    BigEndian::read_u32_into(&hex::decode(val)?, &mut s);
                    for j in 0..L {
                        for k in 0..N {
                            tv.s[j][k] = s[i];
                            i += 1;
                        }
                    }
                },
                "y" => {
                    let mut y = [0; L * N];
                    let mut i = 0;
                    BigEndian::read_u32_into(&hex::decode(val)?, &mut y);
                    for j in 0..L {
                        for k in 0..N {
                            tv.y[j][k] = y[i];
                            i += 1;
                        }
                    }
                },
                "w1" => {
                    let mut w1 = [0; K * N];
                    let mut i = 0;
                    BigEndian::read_u32_into(&hex::decode(val)?, &mut w1);
                    for j in 0..K {
                        for k in 0..N {
                            tv.w1[j][k] = w1[i];
                            i += 1;
                        }
                    }
                },
                "c" => {
                    BigEndian::read_u32_into(&hex::decode(val)?, &mut tv.c);
                },
                _ => panic!()
            }
        }

        testvectors.push(tv);
    }

    Ok(testvectors)
}

#[test]
fn test_vector() {
    for tv in parse_testvectors().unwrap() {
        let mut mat = [PolyVecL::default(); K];
        let mut s = PolyVecL::default();
        let mut y = PolyVecL::default();
        let mut w = PolyVecK::default();
        let mut w1 = PolyVecK::default();
        let mut tmp = PolyVecK::default();
        let mut c = [0; N];

        sign::expand_mat(&mut mat, &tv.seed.0);
        assert!(&mat == &tv.mat);

        for i in 0..L {
            poly::uniform_eta(&mut s[i], &tv.seed.0, i as u8);
        }
        assert!(&s == &tv.s);

        for i in 0..L {
            poly::uniform_gamma1m1(&mut y[i], &tv.seed.0, &tv.seed.1, i as u16);
        }
        assert!(&y == &tv.y);

        y.ntt();
        for i in 0..K {
            polyvec::pointwise_acc_invmontgomery(&mut w[i], &mat[i], &y);
            poly::invntt_montgomery(&mut w[i]);
        }
        w.freeze();
        w.decompose(&mut tmp, &mut w1);
        assert!(&w1 == &tv.w1);

        sign::challenge(&mut c, &tv.seed.1, &w1);
        assert!(&c[..] == &tv.c[..]);
    }
}

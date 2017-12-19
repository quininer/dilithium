#![allow(dead_code)]

use ::params::{ N, L, K };
use ::poly::{ self, Poly };
use ::rounding::{ self, power2round, decompose };
use ::reduce;


macro_rules! polyvec {
    ( $polyvec:ident, $len:expr ) => {
        #[derive(Copy, Clone)]
        pub struct $polyvec(pub [Poly; $len]);

        impl $polyvec {
            pub fn freeze(&mut self) {
                self.0.iter_mut()
                    .for_each(poly::freeze)
            }

            pub fn with_add(&mut self, u: &Self, v: &Self) {
                for i in 0..$len {
                    poly::add(&mut self[i], &u[i], &v[i]);
                }
            }

            pub fn add_assign(&mut self, u: &Self) {
                for i in 0..$len {
                    poly::add_assign(&mut self[i], &u[i]);
                }
            }

            pub fn with_sub(&mut self, u: &Self, v: &Self) {
                for i in 0..$len {
                    poly::sub(&mut self[i], &u[i], &v[i]);
                }
            }

            pub fn neg(&mut self) {
                self.0.iter_mut()
                    .for_each(poly::neg);
            }

            pub fn shift_left(&mut self, k: u32) {
                self.0.iter_mut()
                    .for_each(|p| poly::shift_left(p, k));
            }

            pub fn ntt(&mut self) {
                self.0.iter_mut()
                    .for_each(poly::ntt);
            }

            pub fn invntt_montgomery(&mut self) {
                self.0.iter_mut()
                    .for_each(poly::invntt_montgomery)
            }

            pub fn chknorm(&self, bound: u32) -> bool {
                self.0.iter()
                    .map(|p| poly::chknorm(p, bound))
                    .fold(false, |x, y| x | y)
            }
        }

        impl ::core::ops::Index<usize> for $polyvec {
            type Output = Poly;

            #[inline(always)]
            fn index(&self, i: usize) -> &Self::Output {
                self.0.index(i)
            }
        }

        impl ::core::ops::IndexMut<usize> for $polyvec {
            #[inline(always)]
            fn index_mut(&mut self, i: usize) -> &mut Self::Output {
                self.0.index_mut(i)
            }
        }

        impl ::core::cmp::PartialEq for $polyvec {
            fn eq(&self, other: &Self) -> bool {
                self.0.iter().zip(&other.0)
                    .flat_map(|(x, y)| x.iter().zip(y.iter()))
                    .find(|&(x, y)| x != y)
                    .is_none()
            }
        }

        impl Eq for $polyvec {}

        impl Default for $polyvec {
            fn default() -> Self {
                $polyvec([[0; N]; $len])
            }
        }
    }
}

polyvec!(PolyVecL, L);
polyvec!(PolyVecK, K);

pub fn pointwise_acc_invmontgomery(w: &mut Poly, u: &PolyVecL, v: &PolyVecL) {
    let mut t = [0; N];

    poly::pointwise_invmontgomery(w, &u[0], &v[0]);

    for i in 1..L {
        poly::pointwise_invmontgomery(&mut t, &u[i], &v[i]);
        poly::add_assign(w, &t);
    }

    for i in 0..N {
        w[i] = reduce::reduce32(w[i]);
    }
}

impl PolyVecK {
    pub fn power2round(&self, v0: &mut Self, v1: &mut Self) {
        for i in 0..K {
            for j in 0..N {
                let (x, y) = power2round(self[i][j]);
                v0[i][j] = x;
                v1[i][j] = y;
            }
        }
    }

    pub fn decompose(&self, v0: &mut Self, v1: &mut Self) {
        for i in 0..K {
            for j in 0..N {
                let (x, y) = decompose(self[i][j]);
                v0[i][j] = x;
                v1[i][j] = y;
            }
        }
    }
}

pub fn make_hint(h: &mut PolyVecK, u: &PolyVecK, v: &PolyVecK) -> usize {
    let mut s = 0;
    for i in 0..K {
        for j in 0..N {
            h[i][j] = rounding::make_hint(u[i][j], v[i][j]);
            s += h[i][j] as usize;
        }
    }
    s
}

pub fn use_hint(w: &mut PolyVecK, u: &PolyVecK, h: &PolyVecK) {
    for i in 0..K {
        for j in 0..N {
            w[i][j] = rounding::use_hint(u[i][j], h[i][j]);
        }
    }
}

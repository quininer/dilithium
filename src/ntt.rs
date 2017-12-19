use itertools::Itertools;
use ::params::{ N, Q, MONT, ZETAS, ZETAS_INV };
use ::reduce::{ montgomery_reduce };


pub fn ntt(p: &mut [u32; N]) {
    let mut k = 1;
    for len in (0..8).map(|level| 1 << level).rev() {
        for start in Itertools::step(0..N, 2 * len) {
            let zeta = u64::from(ZETAS[k]);
            k += 1;

            for j in start..(start + len) {
                let t = montgomery_reduce(zeta * u64::from(p[j + len]));
                p[j + len] = p[j] + 2 * Q - t;
                p[j] += t;
            }
        }
    }
}

pub fn invntt_frominvmont(p: &mut [u32; N]) {
    let f = (MONT as u64) * (MONT as u64) % u64::from(Q);
    let f = f * u64::from(Q - 1) % u64::from(Q);
    let f = f * (u64::from(Q - 1) >> 8) % u64::from(Q);

    let mut k = 1;
    for len in (0..8).map(|level| 1 << level) {
        for start in Itertools::step(0..N, 2 * len) {
            let zeta = u64::from(ZETAS_INV[k]);
            k += 1;

            for j in start..(start + len) {
                let t = p[j];
                p[j] += p[j + len];
                p[j + len] = t + 256 * Q - p[j + len];
                p[j + len] = montgomery_reduce(zeta * u64::from(p[j + len]));
            }
        }
    }

    for j in 0..N {
        p[j] = montgomery_reduce(f * u64::from(p[j]));
    }
}

extern crate rand;

use super::*;
use poly::Poly;
use params::{ N, Q };
use self::rand::{ RngCore, thread_rng };


const NTESTS: usize = 10000;

fn poly_naivemul(c: &mut Poly, a: &Poly, b: &Poly) {
    let mut r = [0; 2 * N];

    for i in 0..N {
        for j in 0..N {
            r[i + j] += ((u64::from(a[i]) * u64::from(b[j])) % u64::from(Q)) as u32;
            r[i + j] %= Q;
        }
    }

    for i in N..(2 * N) {
        r[i - N] = r[i - N] + (Q as u32) - r[i];
        r[i - N] %= Q;
    }

    c.copy_from_slice(&r[..N]);
}


#[test]
fn test_mul() {
    let mut rndbuf = [0; 840];
    let (mut c1, mut c2) = ([0; N], [0; N]);
    let (mut a, mut b) = ([0; N], [0; N]);

    for _ in 0..NTESTS {
        thread_rng().fill_bytes(&mut rndbuf);
        poly::uniform(&mut a, &rndbuf);
        thread_rng().fill_bytes(&mut rndbuf);
        poly::uniform(&mut b, &rndbuf);

        poly_naivemul(&mut c1, &a, &b);

        poly::ntt(&mut a);
        poly::ntt(&mut b);
        poly::pointwise_invmontgomery(&mut c2, &a, &b);
        poly::invntt_montgomery(&mut c2);
        poly::csubq(&mut c2);

        assert_eq!(&c2[..], &c1[..]);
    }
}

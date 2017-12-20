use rand::random;
use super::*;
use poly::Poly;
use params::{ N, Q };


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

fn random_poly() -> Poly {
    let mut p = [0; N];
    let mut i = 0;

    while i < N {
        let t = random::<u32>() & 0x7f_ffff;
        if t < Q {
            p[i] = t;
            i += 1;
        }
    }

    p
}


#[test]
fn test_mul() {
    for _ in 0..NTESTS {
        let (mut c1, mut c2) = ([0; N], [0; N]);
        let (mut a, mut b) = (random_poly(), random_poly());

        poly_naivemul(&mut c1, &a, &b);

        poly::ntt(&mut a);
        poly::ntt(&mut b);
        for i in 0..N {
            c2[i] = reduce::montgomery_reduce(u64::from(a[i]) * u64::from(b[i]));
        }
        poly::invntt_montgomery(&mut c2);

        for i in 0..N {
            assert_eq!(reduce::freeze(c2[i]), c1[i]);
        }
    }
}

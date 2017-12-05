use byteorder::{ ByteOrder, LittleEndian };
use ::params::{
    Q, N, ETA, GAMMA1, GAMMA2,
    SEEDBYTES, CRHBYTES
};
use ::reduce;
pub use ::ntt::{ ntt, invntt_frominvmont as invntt_montgomery };


pub type Poly = [u32; N];


pub fn freeze(a: &mut Poly) {
    for i in 0..N {
        a[i] = reduce::freeze(a[i]);
    }
}

pub fn add(c: &mut Poly, a: &Poly, b: &Poly) {
    for i in 0..N {
        c[i] = a[i] + b[i];
    }
}

pub fn add_assign(c: &mut Poly, a: &Poly) {
    for i in 0..N {
        c[i] += a[i];
    }
}

pub fn sub(c: &mut Poly, a: &Poly, b: &Poly) {
    for i in 0..N {
        c[i] = a[i] + 2 * Q -  b[i];
    }
}

pub fn neg(a: &mut Poly) {
    for i in 0..N {
        a[i] = 2 * Q - a[i];
    }
}

pub fn shift_left(a: &mut Poly, k: u32) {
    for i in 0..N {
        a[i] <<= k;
    }
}

pub fn pointwise_invmontgomery(c: &mut Poly, a: &Poly, b: &Poly) {
    for i in 0..N {
        c[i] = reduce::montgomery_reduce(u64::from(a[i]) * u64::from(b[i]));
    }
}

pub fn chknorm(a: &Poly, b: u32) -> bool {
    a.iter()
        .map(|a|{
            let mut t = (Q - 1) / 2 - a;
            t ^= t >> 31;
            (Q - 1) / 2 - t
        })
        .find(|&t| t >= b)
        .is_some()
}

pub fn uniform_eta(a: &mut Poly, seed: &[u8; SEEDBYTES], nonce: u8) {
    use digest::{ Input, ExtendableOutput, XofReader };
    use sha3::Shake256;

    const SHAKE256_RATE: usize = 136;

    fn rej_eta(a: &mut [u32], buf: &[u8]) -> usize {
        let mut ctr = 0;
        let mut pos = 0;

        while ctr < a.len() {
            let (t0, t1) =
                if ETA <= 3 { (u32::from(buf[pos] & 0x07), u32::from(buf[pos] >> 5)) }
                else { (u32::from(buf[pos] & 0x0f), u32::from(buf[pos] >> 5)) };
            pos += 1;

            if t0 <= 2 * ETA {
                a[ctr] = Q + ETA - t0;
                ctr += 1;
            }
            if t1 <= 2 * ETA && ctr < N {
                a[ctr] = Q + ETA - t1;
                ctr += 1;
            }

            if pos >= buf.len() {
                break
            }
        }

        ctr
    }

    let mut outbuf = [0; 2 * SHAKE256_RATE];
    let mut hasher = Shake256::default();
    hasher.process(seed);
    hasher.process(&[nonce]);

    let mut xof = hasher.xof_result();
    xof.read(&mut outbuf);

    let ctr = rej_eta(a, &outbuf);
    if ctr < N {
        xof.read(&mut outbuf[..SHAKE256_RATE]);
        rej_eta(&mut a[ctr..], &outbuf[..SHAKE256_RATE]);
    }
}

pub fn uniform_gamma1m1(a: &mut Poly, seed: &[u8; SEEDBYTES], mu: &[u8; CRHBYTES], nonce: u16) {
    use digest::{ Input, ExtendableOutput, XofReader };
    use sha3::Shake256;

    const SHAKE256_RATE: usize = 136;

    fn rej_gemma1m1(a: &mut [u32], buf: &[u8]) -> usize {
        let mut ctr = 0;
        let mut pos = 0;

        while ctr < a.len() {
            let mut t0 = u32::from(buf[pos]);
            t0 |= u32::from(buf[pos + 1]) << 8;
            t0 |= u32::from(buf[pos + 2]) << 16;
            t0 &= 0xfffff;

            let mut t1 = u32::from(buf[pos + 2]) >> 4;
            t1 |= u32::from(buf[pos + 3]) << 4;
            t1 |= u32::from(buf[pos + 4]) << 12;

            pos += 5;

            if t0 <= 2 * GAMMA1 -2 {
                a[ctr] = Q + GAMMA1 - 1 - t0;
                ctr += 1;
            }
            if t1 <= 2 * GAMMA1 -2 && ctr < a.len() {
                a[ctr] = Q + GAMMA1 - 1 - t1;
                ctr += 1;
            }

            if pos > buf.len() - 5 {
                break
            }
        }

        ctr
    }

    let mut outbuf = [0; 5 * SHAKE256_RATE];
    let mut nonce_bytes = [0; 2];
    LittleEndian::write_u16(&mut nonce_bytes, nonce);

    let mut hasher = Shake256::default();
    hasher.process(seed);
    hasher.process(mu);
    hasher.process(&nonce_bytes);

    let mut xof = hasher.xof_result();
    xof.read(&mut outbuf);

    let ctr = rej_gemma1m1(a, &outbuf);
    if ctr < N {
        xof.read(&mut outbuf[..SHAKE256_RATE]);
        rej_gemma1m1(&mut a[ctr..], &outbuf[..SHAKE256_RATE]);
    }
}

pub fn t1_pack(r: &mut [u8], a: &Poly) {
    for i in 0..(N / 8) {
        r[9*i+0]  = ( a[8*i+0] & 0xFF) as u8;
        r[9*i+1]  = ((a[8*i+0] >> 8) | ((a[8*i+1] & 0x7F) << 1)) as u8;
        r[9*i+2]  = ((a[8*i+1] >> 7) | ((a[8*i+2] & 0x3F) << 2)) as u8;
        r[9*i+3]  = ((a[8*i+2] >> 6) | ((a[8*i+3] & 0x1F) << 3)) as u8;
        r[9*i+4]  = ((a[8*i+3] >> 5) | ((a[8*i+4] & 0x0F) << 4)) as u8;
        r[9*i+5]  = ((a[8*i+4] >> 4) | ((a[8*i+5] & 0x07) << 5)) as u8;
        r[9*i+6]  = ((a[8*i+5] >> 3) | ((a[8*i+6] & 0x03) << 6)) as u8;
        r[9*i+7]  = ((a[8*i+6] >> 2) | ((a[8*i+7] & 0x01) << 7)) as u8;
        r[9*i+8]  = ( a[8*i+7] >> 1) as u8;
    }
}

pub fn t1_unpack(r: &mut Poly, a: &[u8]) {
    for i in 0..(N / 8) {
        r[8*i+0] =  u32::from(a[9*i+0])       | (u32::from(a[9*i+1] & 0x01) << 8);
        r[8*i+1] = (u32::from(a[9*i+1]) >> 1) | (u32::from(a[9*i+2] & 0x03) << 7);
        r[8*i+2] = (u32::from(a[9*i+2]) >> 2) | (u32::from(a[9*i+3] & 0x07) << 6);
        r[8*i+3] = (u32::from(a[9*i+3]) >> 3) | (u32::from(a[9*i+4] & 0x0F) << 5);
        r[8*i+4] = (u32::from(a[9*i+4]) >> 4) | (u32::from(a[9*i+5] & 0x1F) << 4);
        r[8*i+5] = (u32::from(a[9*i+5]) >> 5) | (u32::from(a[9*i+6] & 0x3F) << 3);
        r[8*i+6] = (u32::from(a[9*i+6]) >> 6) | (u32::from(a[9*i+7] & 0x7F) << 2);
        r[8*i+7] = (u32::from(a[9*i+7]) >> 7) | (u32::from(a[9*i+8] & 0xFF) << 1);
    }
}

pub fn w1_pack(r: &mut [u8], a: &Poly) {
    for i in 0..N {
        r[i] = (a[2 * i] | (a[2 * i + 1] << 4)) as u8;
    }
}

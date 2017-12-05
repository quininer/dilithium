use rand::Rng;
use byteorder::{ ByteOrder, LittleEndian };
use ::params::{
    N, K, L, Q,
    BYTES,
    SEEDBYTES, CRHBYTES, POLW1_SIZE_PACKED,
    PK_SIZE_PACKED, SK_SIZE_PACKED, SIG_SIZE_PACKED
};
use ::polyvec::{ self, PolyVecL, PolyVecK };
use ::poly::{ self, Poly };
use ::packing;



/// NOTE panic possible
fn expand_mat(mat: &mut [PolyVecL; K], rho: &[u8; SEEDBYTES]) {
    const SHAKE128_RATE: usize = 168;

    let mut outbuf = [0; 5 * SHAKE128_RATE];

    for i in 0..K {
        for j in 0..L {
            let mut ctr = 0;
            let mut pos = 0;

            shake128!(&mut outbuf; rho, &[(i + (j << 4)) as u8]);

            while ctr < N {
                let val = LittleEndian::read_u24(&outbuf[pos..]) & 0x7fffff;
                pos += 3;

                if val < Q {
                    mat[i][j][ctr] = val;
                    ctr += 1;
                }
            }
        }
    }
}

fn challenge(c: &mut Poly, mu: &[u8; CRHBYTES], w1: &PolyVecK) {
    use digest::{ Input, ExtendableOutput, XofReader };
    use sha3::Shake256;

    const SHAKE256_RATE: usize = 136;

    // zeroed
    c.copy_from_slice(&[0; N]);

    let mut outbuf = [0; SHAKE256_RATE];
    let mut w1pack = [0; K * POLW1_SIZE_PACKED];
    for (i, pack) in w1pack.chunks_mut(POLW1_SIZE_PACKED).enumerate() {
        poly::w1_pack(pack, &w1[i]);
    }

    let mut hasher = Shake256::default();
    hasher.process(mu);
    hasher.process(&w1pack);
    let mut xof = hasher.xof_result();
    xof.read(&mut outbuf);

    let signs = LittleEndian::read_u64(&outbuf);
    let mut pos = 8;
    let mut mask = 1;

    for i in 196..256 {
        let b = loop {
            if pos >= SHAKE256_RATE {
                xof.read(&mut outbuf);
                pos = 0;
            }

            let b = outbuf[pos] as usize;
            pos += 1;
            if b <= i { break b }
        };

        c[i] = c[b];
        c[b] = if signs & mask != 0 { Q - 1 } else { 1 };
        mask <<= 1;
    }
}

pub fn keypair(rng: &mut Rng, pk_bytes: &mut [u8; PK_SIZE_PACKED], sk_bytes: &mut [u8; SK_SIZE_PACKED]) {
    let mut nonce = 0;
    let mut tr = [0; CRHBYTES];
    let mut seedbuf = [0; 3 * SEEDBYTES];
    let mut mat = [PolyVecL::default(); K];
    let mut s1 = PolyVecL::default();
    let (mut s2, mut t, mut t0, mut t1) =
        (PolyVecK::default(), PolyVecK::default(), PolyVecK::default(), PolyVecK::default());

    rng.fill_bytes(&mut seedbuf[..SEEDBYTES]);
    shake256!(&mut seedbuf; &seedbuf[..SEEDBYTES]);
    let rho = array_ref!(seedbuf, 0, SEEDBYTES);
    let rhoprime = array_ref!(seedbuf, SEEDBYTES, SEEDBYTES);
    let key = array_ref!(seedbuf, 2 * SEEDBYTES, SEEDBYTES);

    expand_mat(&mut mat, rho);

    for i in 0..L {
        poly::uniform_eta(&mut s1[i], rhoprime, nonce);
        nonce += 1;
    }
    for i in 0..K {
        poly::uniform_eta(&mut s2[i], rhoprime, nonce);
        nonce += 1;
    }

    let mut s1hat = s1.clone();
    s1hat.ntt();
    for i in 0..K {
        polyvec::pointwise_acc_invmontgomery(&mut t[i], &mat[i], &s1hat);
        poly::invntt_montgomery(&mut t[i])
    }

    t.add_assign(&s2);

    t.freeze();
    t.power2round(&mut t0, &mut t1);
    packing::pk::pack(pk_bytes, &t1, rho);

    shake256!(&mut tr; pk_bytes);
    packing::sk::pack(sk_bytes, rho, key, &tr, &s1, &s2, &t0);
}

pub fn sign(sm: &mut [u8], m: &[u8], sk: &[u8; SK_SIZE_PACKED]) {
    let mut nonce = 0;
    let mut c = [0; N];
    let mut mat = [PolyVecL::default(); K];
    let (mut s1, mut y, mut z) =
        (PolyVecL::default(), PolyVecL::default(), PolyVecL::default());
    let (mut s2, mut t0, mut w, mut w1) =
        (PolyVecK::default(), PolyVecK::default(), PolyVecK::default(), PolyVecK::default());
    let mut tmp = Default::default();
    let (mut rho, mut key, mut mu) = ([0; SEEDBYTES], [0; SEEDBYTES], [0; CRHBYTES]);

    packing::sk::unpack(sk, &mut rho, &mut key, &mut mu, &mut s1, &mut s2, &mut t0);
    shake256!(&mut mu; m, &mu);

    expand_mat(&mut mat, &rho);
    s1.ntt();
    s2.ntt();
    t0.ntt();

    loop {
        for i in 0..L {
            poly::uniform_gamma1m1(&mut y[i], &key, &mu, nonce);
            nonce += 1;
        }

        let mut yhat = y.clone();
        yhat.ntt();
        for i in 0..K {
            polyvec::pointwise_acc_invmontgomery(&mut w[i], &mat[i], &yhat);
            poly::invntt_montgomery(&mut w[i]);
        }

        w.freeze();
        w.decompose(&mut w1, &mut tmp);
        challenge(&mut c, &mu, &w1);

        // TODO
    }
}
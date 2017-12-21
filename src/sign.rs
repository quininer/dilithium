use rand::Rng;
use byteorder::{ ByteOrder, LittleEndian };
use ::params::{
    N, K, L, Q, D, GAMMA1, GAMMA2, BETA, OMEGA,
    SEEDBYTES, CRHBYTES, POLW1_SIZE_PACKED,
    PK_SIZE_PACKED, SK_SIZE_PACKED, SIG_SIZE_PACKED
};
use ::polyvec::{ self, PolyVecL, PolyVecK };
use ::poly::{ self, Poly };
use ::packing;



/// NOTE panic possible
pub(crate) fn expand_mat(mat: &mut [PolyVecL; K], rho: &[u8; SEEDBYTES]) {
    const SHAKE128_RATE: usize = 168;

    let mut outbuf = [0; 5 * SHAKE128_RATE];

    for i in 0..K {
        for j in 0..L {
            let mut ctr = 0;
            let mut pos = 0;

            shake128!(&mut outbuf; rho, &[(i + (j << 4)) as u8]);

            while ctr < N {
                let val = LittleEndian::read_u24(&outbuf[pos..]) & 0x7f_ffff;
                pos += 3;

                if val < Q {
                    mat[i][j][ctr] = val;
                    ctr += 1;
                }
            }
        }
    }
}

pub(crate) fn challenge(c: &mut Poly, mu: &[u8; CRHBYTES], w1: &PolyVecK) {
    use digest::{ Input, ExtendableOutput, XofReader };
    use sha3::Shake256;

    const SHAKE256_RATE: usize = 136;

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

    // Expand 32 bytes of randomness into rho, rhoprime and key
    rng.fill_bytes(&mut seedbuf[..SEEDBYTES]);
    shake256!(&mut seedbuf; &seedbuf[..SEEDBYTES]);
    let rho = array_ref!(seedbuf, 0, SEEDBYTES);
    let rhoprime = array_ref!(seedbuf, SEEDBYTES, SEEDBYTES);
    let key = array_ref!(seedbuf, 2 * SEEDBYTES, SEEDBYTES);

    // Expand matrix
    expand_mat(&mut mat, rho);

    // Sample short vectors s1 and s2
    for i in 0..L {
        poly::uniform_eta(&mut s1[i], rhoprime, nonce);
        nonce += 1;
    }
    for i in 0..K {
        poly::uniform_eta(&mut s2[i], rhoprime, nonce);
        nonce += 1;
    }

    // Matrix-vector multiplication
    let mut s1hat = s1.clone();
    s1hat.ntt();
    for i in 0..K {
        polyvec::pointwise_acc_invmontgomery(&mut t[i], &mat[i], &s1hat);
        poly::invntt_montgomery(&mut t[i])
    }

    // Add noise vector s2
    t.add_assign(&s2);

    // Extract t1 and write public key
    t.freeze();
    t.power2round(&mut t0, &mut t1);
    packing::pk::pack(pk_bytes, &t1, rho);

    // Compute CRH(rho, t1) and write secret key
    shake256!(&mut tr; pk_bytes);
    packing::sk::pack(sk_bytes, rho, key, &tr, &s1, &s2, &t0);


    let mut rho2 = [0; SEEDBYTES];
    let mut key2 = [0; SEEDBYTES];
    let mut mu2 = [0; CRHBYTES];
    let mut s12 = PolyVecL::default();
    let mut s22 = PolyVecK::default();
    let mut t02 = PolyVecK::default();

    packing::sk::unpack(sk_bytes, &mut rho2, &mut key2, &mut mu2, &mut s12, &mut s22, &mut t02);

    assert!(&rho[..] == &rho2[..]);
    assert!(&key[..] == &key2[..]);
    assert!(&mu2[..] == &mu2[..]);
    assert!(s1 == s12);
    assert!(s2 == s22);
    assert!(t0 == t02);
}

pub fn sign(sig: &mut [u8; SIG_SIZE_PACKED], m: &[u8], sk: &[u8; SK_SIZE_PACKED]) {
    let mut nonce = 0;
    let mut mat = [PolyVecL::default(); K];
    let (mut s1, mut y, mut z) =
        (PolyVecL::default(), PolyVecL::default(), PolyVecL::default());
    let (mut s2, mut t0, mut w, mut w1) =
        (PolyVecK::default(), PolyVecK::default(), PolyVecK::default(), PolyVecK::default());
    let (mut h, mut wcs2, mut wcs20, mut ct0, mut tmp) =
        (PolyVecK::default(), PolyVecK::default(), PolyVecK::default(), PolyVecK::default(), PolyVecK::default());
    let (mut rho, mut key, mut mu) = ([0; SEEDBYTES], [0; SEEDBYTES], [0; CRHBYTES]);

    packing::sk::unpack(sk, &mut rho, &mut key, &mut mu, &mut s1, &mut s2, &mut t0);

    // Compute CRH(tr, msg)
    shake256!(&mut mu; &mu, m);

    // Expand matrix and transform vectors
    expand_mat(&mut mat, &rho);
    s1.ntt();
    s2.ntt();
    t0.ntt();

    loop {
        let mut c = [0; N];

        // Sample intermediate vector
        for i in 0..L {
            poly::uniform_gamma1m1(&mut y[i], &key, &mu, nonce);
            nonce += 1;
        }

        // Matrix-vector multiplicatio
        let mut yhat = y.clone();
        yhat.ntt();
        for i in 0..K {
            polyvec::pointwise_acc_invmontgomery(&mut w[i], &mat[i], &yhat);
            poly::invntt_montgomery(&mut w[i]);
        }

        // Decompose w and call the random oracle
        w.freeze();
        w.decompose(&mut tmp, &mut w1);
        challenge(&mut c, &mu, &w1);

        // Compute z, reject if it reveals secret
        let mut chat = c.clone();
        poly::ntt(&mut chat);
        for i in 0..L {
            poly::pointwise_invmontgomery(&mut z[i], &chat, &s1[i]);
            poly::invntt_montgomery(&mut z[i])
        }
        z.add_assign(&y);
        z.freeze();
        if z.chknorm(GAMMA1 - BETA) { continue };

        // Compute w - cs2, reject if w1 can not be computed from it
        for i in 0..K {
            poly::pointwise_invmontgomery(&mut wcs20[i], &chat, &s2[i]);
            poly::invntt_montgomery(&mut wcs20[i]);
        }
        wcs2.with_sub(&w, &wcs20);
        wcs2.freeze();
        wcs2.decompose(&mut wcs20, &mut tmp);
        wcs20.freeze();
        if wcs20.chknorm(GAMMA2 - BETA) { continue };

        if tmp != w1 { continue };

        // Compute hints for w1
        for i in 0..K {
            poly::pointwise_invmontgomery(&mut ct0[i], &chat, &t0[i]);
            poly::invntt_montgomery(&mut ct0[i]);
        }

        ct0.freeze();
        if ct0.chknorm(GAMMA2 - BETA) { continue };

        tmp.with_add(&wcs2, &ct0);
        ct0.neg();
        tmp.freeze();
        let hint = polyvec::make_hint(&mut h, &tmp, &ct0);
        if hint > OMEGA { continue };

        // Write signature
        packing::sign::pack(sig, &z, &h, &c);

        break
    }
}

pub fn verify(m: &[u8], sig: &[u8; SIG_SIZE_PACKED], pk: &[u8; PK_SIZE_PACKED]) -> bool {
    let (mut rho, mut mu) = ([0; SEEDBYTES], [0; CRHBYTES]);
    let (mut c, mut cp) = ([0; N], [0; N]);
    let mut mat = [PolyVecL::default(); K];
    let mut z = PolyVecL::default();
    let (mut t1, mut w1, mut h) = Default::default();
    let (mut tmp1, mut tmp2) = (PolyVecK::default(), PolyVecK::default());

    packing::pk::unpack(pk, &mut t1, &mut rho);
    packing::sign::unpack(sig, &mut z, &mut h, &mut c);
    if z.chknorm(GAMMA1 - BETA) { return false };

    // Compute CRH(CRH(rho, t1), msg)
    shake256!(&mut mu; pk);
    shake256!(&mut mu; &mu, m);

    expand_mat(&mut mat, &rho);

    // Matrix-vector multiplication; compute Az - c2^dt1
    z.ntt();
    for i in 0..K {
        polyvec::pointwise_acc_invmontgomery(&mut tmp1[i], &mat[i], &z);
    }

    let mut chat = c.clone();
    poly::ntt(&mut chat);
    t1.shift_left(D as u32);
    t1.ntt();
    for i in 0..K {
        poly::pointwise_invmontgomery(&mut tmp2[i], &chat, &t1[i]);
    }

    let mut tmp = PolyVecK::default();
    tmp.with_sub(&tmp1, &tmp2);
    tmp.freeze();
    tmp.invntt_montgomery();

    // Reconstruct w1
    tmp.freeze();
    polyvec::use_hint(&mut w1, &tmp, &h);

    // Call random oracle and verify challenge
    challenge(&mut cp, &mu, &w1);
    (0..N)
        .map(|i| c[i] ^ cp[i])
        .fold(0, |sum, next| sum | next)
        .eq(&0)
}

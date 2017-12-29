use ::params::{
    L, K, N, Q, OMEGA,
    SEEDBYTES, CRHBYTES,
    POLETA_SIZE_PACKED, POLT0_SIZE_PACKED, POLT1_SIZE_PACKED, POLZ_SIZE_PACKED,
    PK_SIZE_PACKED, SK_SIZE_PACKED, SIG_SIZE_PACKED
};
use ::poly::{ self, Poly };
use ::polyvec::{ PolyVecL, PolyVecK };


pub mod pk {
    use super::*;

    pub fn pack(pk: &mut [u8; PK_SIZE_PACKED], t1: &PolyVecK, rho: &[u8; SEEDBYTES]) {
        let (rho_bytes, t1_bytes) = mut_array_refs!(pk, SEEDBYTES, POLT1_SIZE_PACKED * K);

        rho_bytes.clone_from(rho);
        for i in 0..K {
            poly::t1_pack(&mut t1_bytes[i * POLT1_SIZE_PACKED..][..POLT1_SIZE_PACKED], &t1[i]);
        }
    }

    pub fn unpack(pk: &[u8; PK_SIZE_PACKED], t1: &mut PolyVecK, rho: &mut [u8; SEEDBYTES]) {
        let (rho_bytes, t1_bytes) = array_refs!(pk, SEEDBYTES, POLT1_SIZE_PACKED * K);

        rho.clone_from(rho_bytes);
        for i in 0..K {
            poly::t1_unpack(&mut t1[i], &t1_bytes[i * POLT1_SIZE_PACKED..][..POLT1_SIZE_PACKED]);
        }
    }
}

pub mod sk {
    use super::*;

    pub fn pack(
        sk: &mut [u8; SK_SIZE_PACKED],
        rho: &[u8; SEEDBYTES],
        key: &[u8; SEEDBYTES],
        tr: &[u8; CRHBYTES],
        s1: &PolyVecL,
        s2: &PolyVecK,
        t0: &PolyVecK
    ) {
        let (rho_bytes, key_bytes, tr_bytes, s1_bytes, s2_bytes, t0_bytes) =
            mut_array_refs!(
                sk,
                SEEDBYTES, SEEDBYTES, CRHBYTES,
                POLETA_SIZE_PACKED * L,
                POLETA_SIZE_PACKED * K,
                POLT0_SIZE_PACKED * K
            );

        rho_bytes.clone_from(rho);
        key_bytes.clone_from(key);
        tr_bytes.clone_from(tr);

        for i in 0..L {
            poly::eta_pack(&mut s1_bytes[i * POLETA_SIZE_PACKED..][..POLETA_SIZE_PACKED], &s1[i]);
        }
        for i in 0..K {
            poly::eta_pack(&mut s2_bytes[i * POLETA_SIZE_PACKED..][..POLETA_SIZE_PACKED], &s2[i]);
        }
        for i in 0..K {
            poly::t0_pack(&mut t0_bytes[i * POLT0_SIZE_PACKED..][..POLT0_SIZE_PACKED], &t0[i]);
        }
    }

    pub fn unpack(
        sk: &[u8; SK_SIZE_PACKED],
        rho: &mut [u8; SEEDBYTES],
        key: &mut [u8; SEEDBYTES],
        tr: &mut [u8; CRHBYTES],
        s1: &mut PolyVecL,
        s2: &mut PolyVecK,
        t0: &mut PolyVecK
   ) {
        let (rho_bytes, key_bytes, tr_bytes, s1_bytes, s2_bytes, t0_bytes) =
            array_refs!(
                sk,
                SEEDBYTES, SEEDBYTES, CRHBYTES,
                POLETA_SIZE_PACKED * L,
                POLETA_SIZE_PACKED * K,
                POLT0_SIZE_PACKED * K
            );

        rho.clone_from(rho_bytes);
        key.clone_from(key_bytes);
        tr.clone_from(tr_bytes);

        for i in 0..L {
            poly::eta_unpack(&mut s1[i], &s1_bytes[i * POLETA_SIZE_PACKED..][..POLETA_SIZE_PACKED]);
        }
        for i in 0..K {
            poly::eta_unpack(&mut s2[i], &s2_bytes[i * POLETA_SIZE_PACKED..][..POLETA_SIZE_PACKED]);
        }
        for i in 0..K {
            poly::t0_unpack(&mut t0[i], &t0_bytes[i * POLT0_SIZE_PACKED..][..POLT0_SIZE_PACKED]);
        }
    }
}

pub mod sign {
    use super::*;

    pub fn pack(sign: &mut [u8; SIG_SIZE_PACKED], z: &PolyVecL, h: &PolyVecK,c: &Poly) {
        let (z_bytes, h_bytes, c_bytes) =
            mut_array_refs!(
                sign,
                POLZ_SIZE_PACKED * L,
                OMEGA + K,
                N / 8 + 8
            );

        for i in 0..L {
            poly::z_pack(&mut z_bytes[i * POLZ_SIZE_PACKED..][..POLZ_SIZE_PACKED], &z[i]);
        }

        let mut k = 0;
        for i in 0..K {
            for j in 0..N {
                if h[i][j] == 1 {
                    h_bytes[k] = j as u8;
                    k += 1;
                }
            }
            h_bytes[OMEGA + i] = k as u8;
        }

        let mut signs: u64 = 0;
        let mut mask = 1;
        for i in 0..(N / 8) {
            for j in 0..8 {
                if c[8 * i + j] != 0 {
                    c_bytes[i] |= 1 << j;
                    if c[8 * i + j] == Q - 1 {
                        signs |= mask;
                    }
                    mask <<= 1;
                }
            }
        }
        for i in 0..8 {
            c_bytes[N / 8..][i] = (signs >> (8 * i)) as u8;
        }
    }

    pub fn unpack(sign: &[u8; SIG_SIZE_PACKED], z: &mut PolyVecL, h: &mut PolyVecK,c: &mut Poly) {
        let (z_bytes, h_bytes, c_bytes) =
            array_refs!(
                sign,
                POLZ_SIZE_PACKED * L,
                OMEGA + K,
                N / 8 + 8
            );
        for i in 0..L {
            poly::z_unpack(&mut z[i], &z_bytes[i * POLZ_SIZE_PACKED..][..POLZ_SIZE_PACKED]);
        }

        let mut k = 0;
        for i in 0..K {
            for j in k..(h_bytes[OMEGA + i] as usize) {
                h[i][h_bytes[j] as usize] = 1;
            }
            k = h_bytes[OMEGA + i] as usize;
        }

        let signs = (0..8)
            .map(|i| u64::from(c_bytes[N / 8 + i]) << (8 * i))
            .fold(0, |sum, next| sum | next);
        let mut mask = 1;
        for i in 0..(N / 8) {
            for j in 0..8 {
                if (c_bytes[i] >> j) & 0x01 != 0 {
                    c[8 * i + j] =
                        if (signs & mask) != 0 { Q - 1 }
                        else { 1 };
                    mask <<= 1;
                }
            }
        }
    }
}

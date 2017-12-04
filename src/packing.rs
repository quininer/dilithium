use ::params::{
    SEEDBYTES, CRHBYTES,
    PK_SIZE_PACKED, SK_SIZE_PACKED, SIG_SIZE_PACKED
};
use ::poly::Poly;
use ::polyvec::{ PolyVecL, PolyVecK };




pub mod pk {
    use super::*;

    pub fn pack(pk: &mut [u8; PK_SIZE_PACKED], t1: &PolyVecK, rho: &[u8; SEEDBYTES]) {
        unimplemented!()
    }

    pub fn unpack(pk: &[u8; PK_SIZE_PACKED], t1: &mut PolyVecK, rho: &mut [u8; SEEDBYTES]) {
        unimplemented!()
    }
}

pub mod sk {
    use super::*;

    pub fn pack(
        sk: &mut [u8; SK_SIZE_PACKED],
        s1: &PolyVecL,
        s2: &PolyVecK,
        t0: &PolyVecK,
        rho: &[u8; SEEDBYTES],
        key: &[u8; SEEDBYTES],
        tr: &[u8; CRHBYTES]
    ) {
        unimplemented!()
    }

    pub fn unpack(
        sk: &[u8; SK_SIZE_PACKED],
        s1: &mut PolyVecL,
        s2: &mut PolyVecK,
        t0: &mut PolyVecK,
        rho: &mut [u8; SEEDBYTES],
        key: &mut [u8; SEEDBYTES],
        tr: &mut [u8; CRHBYTES]
    ) {
        unimplemented!()
    }
}

pub mod sign {
    use super::*;

    pub fn pack(
        sign: &mut [u8; SIG_SIZE_PACKED],
        z: &PolyVecL,
        h: &PolyVecK,
        c: &Poly
    ) {
        unimplemented!()
    }

    pub fn unpack(
        sign: &[u8; SIG_SIZE_PACKED],
        z: &mut PolyVecL,
        h: &mut PolyVecK,
        c: &mut Poly
    ) {
        unimplemented!()
    }
}

use ::params::{
    SEEDBYTES, CRHBYTES,
    PK_SIZE_PACKED, SK_SIZE_PACKED, SIG_SIZE_PACKED
};
use ::poly::Poly;
use ::polyvec::{ PolyVecL, PolyVecK };


pub struct PublicKey {
    t1: PolyVecK,
    rho: [u8; SEEDBYTES]
}

pub struct SecretKey {
    s1: PolyVecL,
    s2: PolyVecK,
    t0: PolyVecK,
    rho: [u8; SEEDBYTES],
    key: [u8; SEEDBYTES],
    tr: [u8; CRHBYTES]
}

pub struct Signature {
    z: PolyVecL,
    h: PolyVecK,
    c: Poly
}


impl PublicKey {
    pub fn pack(&self, pk: &mut [u8; PK_SIZE_PACKED]) {
        unimplemented!()
    }

    pub fn unpack(&mut self, pk: &[u8; PK_SIZE_PACKED]) {
        unimplemented!()
    }
}

impl SecretKey {
    pub fn pack(&self, sk: &mut [u8; SK_SIZE_PACKED]) {
        unimplemented!()
    }

    pub fn unpack(&mut self, sk: &[u8; SK_SIZE_PACKED]) {
        unimplemented!()
    }
}

impl Signature {
    pub fn pack(&self, sign: &mut [u8; SIG_SIZE_PACKED]) {
        unimplemented!()
    }

    pub fn unpack(&mut self, sign: &[u8; SIG_SIZE_PACKED]) {
        unimplemented!()
    }
}

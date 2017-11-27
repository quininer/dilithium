pub const SEEDBYTES    : usize = 32;
pub const CRHBYTES     : usize = 48;
pub const N            : usize = 256;
pub const Q            : usize = 8380417;
pub const QBITS        : usize = 23;
pub const ROOT_OF_UNITY: usize = 1753;
pub const D            : usize = 14;
pub const GAMMA1       : usize = (Q - 1) / 16;
pub const GAMMA2       : usize = GAMMA1 / 2;
pub const ALPHA        : usize = 2 * GAMMA2;

#[cfg(feature = "mode0")]
mod mode {
    pub const K       : usize = 3;
    pub const L       : usize = 2;
    pub const ETA     : usize = 7;
    pub const SETABITS: usize = 4;
    pub const BETA    : usize =375;
    pub const OMEGA   : usize = 64;
}


#[cfg(feature = "mode1")]
mod mode {
    pub const K       : usize = 4;
    pub const L       : usize = 3;
    pub const ETA     : usize = 6;
    pub const SETABITS: usize = 4;
    pub const BETA    : usize = 325;
    pub const OMEGA   : usize = 80;
}

#[cfg(feature = "mode2")]
mod mode {
    pub const K       : usize = 5;
    pub const L       : usize = 4;
    pub const ETA     : usize = 5;
    pub const SETABITS: usize = 4;
    pub const BETA    : usize = 275;
    pub const OMEGA   : usize = 96;
}

#[cfg(feature = "mode3")]
mod mode {
    pub const K       : usize = 6;
    pub const L       : usize = 5;
    pub const ETA     : usize = 3;
    pub const SETABITS: usize = 3;
    pub const BETA    : usize = 175;
    pub const OMEGA   : usize = 120;
}

pub use self::mode::*;


pub const POL_SIZE_PACKED   : usize = (N * QBITS) / 8;
pub const POLT1_SIZE_PACKED : usize = (N * (QBITS - D)) / 8;
pub const POLT0_SIZE_PACKED : usize = (N * D) / 8;
pub const POLETA_SIZE_PACKED: usize = (N * SETABITS) / 8;
pub const POLZ_SIZE_PACKED  : usize = (N * (QBITS - 3)) / 8;
pub const POLW1_SIZE_PACKED : usize = (N * 4) / 8;

pub const POLVECK_SIZE_PACKED: usize = K * POL_SIZE_PACKED;
pub const POLVECL_SIZE_PACKED: usize = L * POL_SIZE_PACKED;
pub const PK_SIZE_PACKED     : usize = SEEDBYTES + K * POLT1_SIZE_PACKED;
pub const SK_SIZE_PACKED     : usize = 2 * SEEDBYTES + (L + K) * POLETA_SIZE_PACKED + CRHBYTES + K * POLT0_SIZE_PACKED;
pub const SIG_SIZE_PACKED    : usize = L * POLZ_SIZE_PACKED + (OMEGA + K) + (N / 8 + 8);


pub const MONT: usize = 4193792;
pub const QINV: usize = 4236238847;

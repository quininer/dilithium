extern crate core;

#[macro_use] extern crate arrayref;
extern crate rand;
extern crate byteorder;
extern crate itertools;
extern crate digest;
extern crate sha3;

#[macro_use] mod utils;
pub mod params;
pub mod reduce;
pub mod rounding;
pub mod ntt;
pub mod poly;
pub mod polyvec;
pub mod packing;
pub mod sign;

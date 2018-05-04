extern crate core;

#[macro_use] extern crate arrayref;
extern crate rand_core;
extern crate byteorder;
extern crate itertools;
extern crate digest;
extern crate sha3;

#[macro_use] mod utils;
mod reduce;
mod rounding;
mod ntt;
mod poly;
mod polyvec;
mod packing;
pub mod params;
pub mod sign;

#[cfg(test)] mod test_mul;
#[cfg(test)] mod test_vectors;

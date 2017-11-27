use ::params::{ Q, QINV };


pub fn montgomery_reduce(a: u64) -> u32 {
    let mut t = a * QINV as u64;
    t &= (1 << 32) - 1;
    t *= Q as u64;
    t += a;
    (t >> 32) as u32
}

pub fn reduce32(mut a: u32) -> u32 {
    let mut t = a & 0x7FFFFF;
    a >>= 23;
    t += (a << 13) - a;
    t
}

pub fn freeze(a: u32) -> u32 {
    let mut a = reduce32(a);
    a = a.wrapping_sub(Q as u32);
    let c = ((a as i32) >> 31) & Q as i32;
    a + (c as u32)
}

use ::params::{ Q, QINV };


pub fn montgomery_reduce(a: u64) -> u32 {
    let mut t = a.wrapping_mul(QINV as u64);
    t &= (1 << 32) - 1;
    t *= u64::from(Q);
    t += a;
    (t >> 32) as u32
}

pub fn reduce32(mut a: u32) -> u32 {
    let mut t = a & 0x7f_ffff;
    a >>= 23;
    t += (a << 13) - a;
    t
}

pub fn csubq(mut a: u32) -> u32 {
    a = a.wrapping_sub(Q as u32);
    let c = ((a as i32) >> 31) & Q as i32;
    a.wrapping_add(c as u32)
}

pub fn freeze(a: u32) -> u32 {
    let a = reduce32(a);
    let a = csubq(a);
    a
}

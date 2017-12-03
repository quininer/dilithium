use ::params::{ Q, D, ALPHA };
use ::reduce::freeze;


pub fn power2round(a: u32) -> (u32, u32) {
    let d = D as u32;

    let mut t = a & ((1 << d) - 1);
    t -= (1 << (d - 1)) + 1;
    t += (t >> 31) & (1 << d);
    t -= (1 << (d - 1)) - 1;
    (Q + t, (a - t) >> d)
}

pub fn decompose(mut a: u32) -> (u32, u32) {
    let mut t = a & 0x7FFFF;
    t += (a >> 19) << 9;
    t -= ALPHA / 2 + 1;
    t += (t >> 31) & ALPHA;
    t -= ALPHA / 2 - 1;
    a -= t;

    let mut u = a - 1;
    u >>= 31;
    a = (a >> 19) + 1;
    a -= u & 1;

    (Q + t - (a >> 4), a & 0xf)
}

pub fn make_hint(a: u32, b: u32) -> u32 {
    let (_, x) = decompose(a);
    let (_, y) = decompose(freeze(a + b));
    if x != y { 1 } else { 0 }
}

pub fn use_hint(a: u32, hint: u32) -> u32 {
    let (a0, a1) = decompose(a);

    if hint == 0 {
        a1
    } else if a0 > Q {
        if a1 == (Q - 1) / ALPHA - 1 {
            0
        } else {
            a1 + 1
        }
    } else {
        if a1 == 0 {
            (Q - 1) / ALPHA - 1
        } else {
            a1 - 1
        }
    }
}

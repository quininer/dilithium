extern crate rand;
extern crate dilithium;

use rand::{ Rng, thread_rng };
use dilithium::params::*;
use dilithium::sign::{ keypair, sign, verify };


#[test]
fn test_sign() {
    let mut rng = thread_rng();
    let mut message = [0; 59];
    let (mut pk, mut sk) = ([0; PK_SIZE_PACKED], [0; SK_SIZE_PACKED]);
    let mut sig = [0; SIG_SIZE_PACKED];

    rng.fill_bytes(&mut message);

    keypair(&mut rng, &mut pk, &mut sk);
    sign(&mut sig, &message, &sk);
    assert!(verify(&message, &sig, &pk));
}

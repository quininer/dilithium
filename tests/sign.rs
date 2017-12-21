extern crate rand;
extern crate dilithium;

use rand::{ Rng, thread_rng };
use dilithium::params::*;
use dilithium::sign::{ keypair, sign, verify };


#[test]
fn test_sign() {
    for _ in 0..5000 {
        let mut rng = thread_rng();
        let mut message = [0; 59];
        let (mut pk, mut sk) = ([0; PUBLICKEYBYTES], [0; SECRETKEYBYTES]);
        let mut sig = [0; BYTES];
        rng.fill_bytes(&mut message);

        keypair(&mut rng, &mut pk, &mut sk);
        sign(&mut sig, &message, &sk);

        assert!(verify(&message, &sig, &pk));

        message[2] ^= 42;
        assert!(!verify(&message, &sig, &pk));
    }
}

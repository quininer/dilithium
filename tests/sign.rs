extern crate rand;
extern crate dilithium;

use rand::{ RngCore, FromEntropy, ChaChaRng };
use dilithium::params::*;
use dilithium::sign::{ keypair, sign, verify };


#[test]
fn test_sign() {
    for _ in 0..500 {
        let mut rng = ChaChaRng::from_entropy();
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

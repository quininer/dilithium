macro_rules! shake128 {
    ( $output:expr ; $( $input:expr ),* ) => {
        let mut hasher = ::sha3::Shake128::default();
        $(
            ::digest::Input::process(&mut hasher, $input);
        )*
        let mut reader = ::digest::ExtendableOutput::xof_result(hasher);
        ::digest::XofReader::read(&mut reader, $output);
    }
}

macro_rules! shake256 {
    ( $output:expr ; $( $input:expr ),* ) => {
        let mut hasher = ::sha3::Shake256::default();
        $(
            ::digest::Input::process(&mut hasher, $input);
        )*
        let mut reader = ::digest::ExtendableOutput::xof_result(hasher);
        ::digest::XofReader::read(&mut reader, $output);
    }
}

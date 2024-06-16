#![feature(portable_simd)]
// Listening to this warning would break most functions
#![allow(clippy::needless_range_loop)]
// This causes issues w/ the self as parameter of methods
#![allow(clippy::only_used_in_recursion)]
#![allow(clippy::new_without_default)]

pub mod algos;

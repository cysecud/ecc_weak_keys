use std::fmt::Debug;

use num_traits::{One, Zero};

#[macro_use]
pub mod ellinit;
pub mod affine;
pub mod projective;
pub mod nistp192;
pub mod nistp224;
pub mod nistp256k1;
pub mod nistp256r1;
pub mod nistsm256;
pub mod nistp384;
pub mod nistp521;

pub trait Scalar<'a> {
    const SIZE:usize;
    const ORDER:&'a str;
}
pub trait FieldElement<'a,T:Scalar<'a>>{
    const PRIME:&'a str;
    const PRIME_ROOT:&'a str;
    const ZERO:&'a str ="0";
    const ONE:&'a str="1";
    fn zero()->Self;
    fn one()->Self;
    fn prime()->T;
    fn new(num:T)->Self; 
    fn double(&mut self)->Self;
    fn multiple(&self,n:T)->Self;
    fn square(&mut self)->Self;
    fn power(&mut self, n:T)->Self;
    fn inverse(&mut self)->Self;
    fn rand_mod()->Self;
}
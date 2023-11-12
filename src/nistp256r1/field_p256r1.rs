use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};

use crate::FieldElement;
use super::scalar256::U256;

#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct FieldP256r1{
    pub num:U256,
}
impl Add for FieldP256r1{
    type Output = Self;
    fn add(self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        FieldP256r1{num:sum.mod_fast(Self::prime())}}
}  
impl Mul for FieldP256r1{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        FieldP256r1{num:mul.mod_fast(Self::prime())}}
}
impl Neg for FieldP256r1{
    type Output = Self;
    fn neg(self) -> Self::Output {
        FieldP256r1 {num:Self::prime()-self.num}
    }
}
impl Sub for FieldP256r1{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
            FieldP256r1 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
    }
}
impl <'a>FieldElement<'a,U256> for FieldP256r1{
    const PRIME:&'a str ="115792089237316195423570985008687907853269984665640564039457584007908834671663" ;
    fn zero()->Self {
        Self { num: U256::zero() }
    }
    fn one()->Self {
        Self { num: U256::one() }
    }
    fn prime()->U256 {
        U256::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U256)->FieldP256r1{
          let prime = U256::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            FieldP256r1 {num:num}
        }
    }
    fn double(&mut self)->FieldP256r1 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U256)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U256::from(2),Self::prime());
        FieldP256r1::new(sqr)
    }
    fn power(&mut self,n: U256) -> Self {
        /*Montgomery ladder */
        
        let bin=U256::to_binary(n);
        let mut p0=FieldP256r1::new(U256::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U256::zero() {
                p1=p1*p0;
                p0=p0.square();
            }
            else {
                p0=p0*p1;
                p1=p1.square();
                }
        }
        p0
        
    }
    fn inverse (&mut self)-> FieldP256r1 {
        FieldP256r1::power(self, Self::prime()-U256::from(2))
}
fn rand_mod()-> FieldP256r1{
    let prime=U256::from_dec_str(Self::PRIME).expect("error");
    let num=U256::random().mod_fast(prime);
    FieldP256r1::new(num)
}
    }

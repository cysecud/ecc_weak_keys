use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};

use crate::FieldElement;
use super::scalar521::U521;

#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitP521{
    pub num:U521,
}
impl Add for ImplicitP521{
    type Output = Self;
    fn add(self:Self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        ImplicitP521{num:sum.mod_fast(Self::prime())}}

}  
impl Mul for ImplicitP521{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        ImplicitP521{num:mul.mod_fast(Self::prime())}}
}
impl Neg for ImplicitP521{
    type Output = Self;
    fn neg(self) -> Self::Output {
        ImplicitP521 {num:Self::prime()-self.num}
    }
}
impl Sub for ImplicitP521{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        
            ImplicitP521 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
        }  
}
impl <'a>FieldElement<'a,U521> for ImplicitP521{
    const PRIME:&'a str ="6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449";
    const PRIME_ROOT:&'a str = "3";
    fn zero()->Self {
        Self { num: U521::zero() }
    }
    fn one()->Self {
        Self { num: U521::one() }
    }
    fn prime()->U521 {
        U521::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U521)->ImplicitP521{
          let prime = U521::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            ImplicitP521 {num:num}
        }
    }
    fn double(&mut self)->ImplicitP521 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U521)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U521::from(2),Self::prime());
        ImplicitP521::new(sqr)
    }
    fn power(&mut self,n: U521) -> Self {
        /*Montgomery ladder */
        
        let bin=U521::to_binary(n);
        let mut p0=ImplicitP521::new(U521::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U521::zero() {
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
    fn inverse (&mut self)-> ImplicitP521 {
        ImplicitP521::power(self, Self::prime()-U521::from(2))
}
fn rand_mod()-> ImplicitP521{
    let prime=U521::from_dec_str(Self::PRIME).expect("error");
    let num=U521::random().mod_fast(prime);
    ImplicitP521::new(num)
}
    }

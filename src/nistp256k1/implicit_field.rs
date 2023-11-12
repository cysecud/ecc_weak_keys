use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};
use super::FieldElement;
use super::scalar256::U256;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitP256k1{
    pub num:U256,
}
impl Add for ImplicitP256k1{
    type Output = Self;
    fn add(self:Self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        ImplicitP256k1{num:sum.mod_fast(Self::prime())}}

}  
impl Mul for ImplicitP256k1{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        ImplicitP256k1{num:mul.mod_fast(Self::prime())}}
}
impl Neg for ImplicitP256k1{
    type Output = Self;
    fn neg(self) -> Self::Output {
        ImplicitP256k1 {num:Self::prime()-self.num}
    }
}
impl Sub for ImplicitP256k1{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        
            ImplicitP256k1 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
        }  
}
impl <'a>FieldElement<'a,U256> for ImplicitP256k1{
    const PRIME:&'a str ="26959946667150639794667015087019625940457807714424391721682722368061";
    fn zero()->Self {
        Self { num: U256::zero() }
    }
    fn one()->Self {
        Self { num: U256::one() }
    }
    fn prime()->U256 {
        U256::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U256)->ImplicitP256k1{
          let prime = U256::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            ImplicitP256k1 {num:num}
        }
    }
    fn double(&mut self)->ImplicitP256k1 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U256)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U256::from(2),Self::prime());
        ImplicitP256k1::new(sqr)
    }
    fn power(&mut self,n: U256) -> Self {
        /*Montgomery ladder */
        
        let bin=U256::to_binary(n);
        let mut p0=ImplicitP256k1::new(U256::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U256::zero() {
                p1=p1*p0;
                p0=p0.square();
                
                /* println!("p0 is {:?}",p0);
                println!("p1 is {:?}",p1); */
            }
            else {
                p0=p0*p1;
                p1=p1.square();
              /*   println!("p0 is {:?}",p0);
                println!("p1 is {:?}",p1); */
                }
        }
        p0
        
    }
    fn inverse (&mut self)-> ImplicitP256k1 {
        ImplicitP256k1::power(self, Self::prime()-U256::from(2))
}
fn rand_mod()-> ImplicitP256k1{
    let prime=U256::from_dec_str(Self::PRIME).expect("error");
    let num=U256::random().mod_fast(prime);
    ImplicitP256k1::new(num)
}
    }

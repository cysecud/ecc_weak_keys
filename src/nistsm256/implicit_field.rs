use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};
use crate::FieldElement;
use super::scalar256::U256;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitSM256{
    pub num:U256,
}
impl Add for ImplicitSM256{
    type Output = Self;
    fn add(self:Self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        ImplicitSM256{num:sum.mod_fast(Self::prime())}}

}  
impl Mul for ImplicitSM256{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        ImplicitSM256{num:mul.mod_fast(Self::prime())}}
}
impl Neg for ImplicitSM256{
    type Output = Self;
    fn neg(self) -> Self::Output {
        ImplicitSM256 {num:Self::prime()-self.num}
    }
}
impl Sub for ImplicitSM256{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        
            ImplicitSM256 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
        }  
}
impl <'a>FieldElement<'a,U256> for ImplicitSM256{
    const PRIME:&'a str ="115792089210356248762697446949407573529996955224135760342422259061068512044369";
    fn zero()->Self {
        Self { num: U256::zero() }
    }
    fn one()->Self {
        Self { num: U256::one() }
    }
    fn prime()->U256 {
        U256::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U256)->ImplicitSM256{
          let prime = U256::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            ImplicitSM256 {num:num}
        }
    }
    fn double(&mut self)->ImplicitSM256 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U256)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U256::from(2),Self::prime());
        ImplicitSM256::new(sqr)
    }
    fn power(&mut self,n: U256) -> Self {
        /*Montgomery ladder */
        
        let bin=U256::to_binary(n);
        let mut p0=ImplicitSM256::new(U256::one());
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
    fn inverse (&mut self)-> ImplicitSM256 {
        ImplicitSM256::power(self, Self::prime()-U256::from(2))
}
fn rand_mod()-> ImplicitSM256{
    let prime=U256::from_dec_str(Self::PRIME).expect("error");
    let num=U256::random().mod_fast(prime);
    ImplicitSM256::new(num)
}
    }

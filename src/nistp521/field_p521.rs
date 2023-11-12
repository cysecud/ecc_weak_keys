use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};

use crate::FieldElement;
use super::scalar521::U521;

#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct FieldP521{
    pub num:U521,
}
impl Add for FieldP521{
    type Output = Self;
    fn add(self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        FieldP521{num:sum.mod_fast(Self::prime())}}
}  
impl Mul for FieldP521{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        FieldP521{num:mul.mod_fast(Self::prime())}}
}
impl Neg for FieldP521{
    type Output = Self;
    fn neg(self) -> Self::Output {
        FieldP521 {num:Self::prime()-self.num}
    }
}
impl Sub for FieldP521{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
            FieldP521 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
    }
}
impl <'a>FieldElement<'a,U521> for FieldP521{
    const PRIME:&'a str ="6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151" ;
    fn zero()->Self {
        Self { num: U521::zero() }
    }
    fn one()->Self {
        Self { num: U521::one() }
    }
    fn prime()->U521 {
        U521::from_dec_str(Self::PRIME).expect("error in Prime Field 521")
    }
    fn new(num:U521)->FieldP521{
          let prime = U521::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            FieldP521 {num:num}
        }
    }
    fn double(&mut self)->FieldP521 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U521)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let  sqr=self.num.pow_mod(U521::from(2),Self::prime());
        FieldP521::new(sqr)
    }
    fn power(&mut self,n: U521) -> Self {
        /*Montgomery ladder */
        
        let bin=U521::to_binary(n);
        let mut p0=FieldP521::new(U521::one());
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
    fn inverse (&mut self)-> FieldP521 {
        FieldP521::power(self, Self::prime()-U521::from(2))
}
fn rand_mod()-> FieldP521{
    let prime=FieldP521::prime();
    let num=U521::random().mod_fast(prime);
    FieldP521::new(num)
}
    }

use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};


use super::FieldElement;
use super::scalar384::U384;

#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct FieldP384{
    pub num:U384,
}
impl Add for FieldP384{
    type Output = Self;
    fn add(self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        FieldP384{num:sum.mod_fast(Self::prime())}}
}  
impl Mul for FieldP384{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        FieldP384{num:mul.mod_fast(Self::prime())}}
}
impl Neg for FieldP384{
    type Output = Self;
    fn neg(self) -> Self::Output {
        FieldP384 {num:Self::prime()-self.num}
    }
}
impl Sub for FieldP384{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
            FieldP384 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
    }
}
impl <'a>FieldElement<'a,U384> for FieldP384{
    const PRIME:&'a str ="39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319" ;
    const PRIME_ROOT:&'a str = "19";
    fn zero()->Self {
        FieldP384 { num: U384::zero() }
    }
    fn one()->Self {
        FieldP384 { num: U384::one() }
    }
    fn prime()->U384 {
        U384::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U384)->FieldP384{
          let prime = U384::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            FieldP384 {num:num}
        }
    }
    fn double(&mut self)->FieldP384 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U384::from(2),Self::prime());
        FieldP384::new(sqr)
    }
    fn multiple(&self,n:U384)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn power(&mut self,n: U384) -> Self {
        /*Montgomery ladder */
        
        let bin=U384::to_binary(n);
        let mut p0=FieldP384::new(U384::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U384::zero() {
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
    fn inverse (&mut self)-> FieldP384 {
        FieldP384::power(self, Self::prime()-U384::from(2))
}
fn rand_mod()-> FieldP384{
    let prime=U384::from_dec_str(Self::PRIME).expect("error");
    let num=U384::random().mod_fast(prime);
    FieldP384::new(num)
}
    }

use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};
use super::FieldElement;
use super::scalar224::U224;

#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct FieldP224{
    pub num:U224,
}
impl Add for FieldP224{
    type Output = Self;
    fn add(self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        FieldP224{num:sum.mod_fast(Self::prime())}}
}  
impl Mul for FieldP224{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        FieldP224{num:mul.mod_fast(Self::prime())}}
}
impl Neg for FieldP224{
    type Output = Self;
    fn neg(self) -> Self::Output {
        FieldP224 {num:Self::prime()-self.num}
    }
}
impl Sub for FieldP224{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
            FieldP224 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
    }
}
impl <'a>FieldElement<'a,U224> for FieldP224{
    const PRIME:&'a str ="26959946667150639794667015087019630673557916260026308143510066298881" ;
    fn zero()->Self {
        Self { num: U224::zero() }
    }
    fn one()->Self {
        Self { num: U224::one() }
    }
    fn prime()->U224 {
        U224::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U224)->FieldP224{
          let prime = U224::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            FieldP224 {num:num}
        }
    }
    fn double(&mut self)->FieldP224 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U224)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U224::from(2),Self::prime());
        FieldP224::new(sqr)
    }
    fn power(&mut self,n: U224) -> Self {
        /*Montgomery ladder */
        
        let bin=U224::to_binary(n);
        let mut p0=FieldP224::new(U224::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U224::zero() {
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
    fn inverse (&mut self)-> FieldP224 {
        FieldP224::power(self, Self::prime()-U224::from(2))
}
fn rand_mod()-> FieldP224{
    let prime=U224::from_dec_str(Self::PRIME).expect("error");
    let num=U224::random().mod_fast(prime);
    FieldP224::new(num)
}
    }

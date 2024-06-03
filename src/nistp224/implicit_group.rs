use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};
use super::FieldElement;
use super::scalar224::U224;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitP224{
    pub num:U224,
}
impl ImplicitP224{
    pub const DIVISOR_32:[&str; 1]=["533642580"];
    /*DIVISOR_64, DIVISOR_128, DIVISOR_160 are the same as DIVISOR_32.
    They are computed by using Pari_gp*/
}

impl Add for ImplicitP224{
    type Output = Self;
    fn add(self:Self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        ImplicitP224{num:sum.mod_fast(Self::prime())}}

}  
impl Mul for ImplicitP224{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        ImplicitP224{num:mul.mod_fast(Self::prime())}}
}
impl Neg for ImplicitP224{
    type Output = Self;
    fn neg(self) -> Self::Output {
        ImplicitP224 {num:Self::prime()-self.num}
    }
}
impl Sub for ImplicitP224{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        
            ImplicitP224 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
        }  
}
impl <'a>FieldElement<'a,U224> for ImplicitP224{
    const PRIME:&'a str ="26959946667150639794667015087019625940457807714424391721682722368061";
    const PRIME_ROOT:&'a str = "2";
    fn zero()->Self {
        Self { num: U224::zero() }
    }
    fn one()->Self {
        Self { num: U224::one() }
    }
    fn prime()->U224 {
        U224::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U224)->ImplicitP224{
          let prime = U224::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            ImplicitP224 {num:num}
        }
    }
    fn double(&mut self)->ImplicitP224 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U224)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U224::from(2),Self::prime());
        ImplicitP224::new(sqr)
    }
    fn power(&mut self,n: U224) -> Self {
        /*Montgomery ladder */
        
        let bin=U224::to_binary(n);
        let mut p0=ImplicitP224::new(U224::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U224::zero() {
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
    fn inverse (&mut self)-> ImplicitP224 {
        ImplicitP224::power(self, Self::prime()-U224::from(2))
}
fn rand_mod()-> ImplicitP224{
    let prime=U224::from_dec_str(Self::PRIME).expect("error");
    let num=U224::random().mod_fast(prime);
    ImplicitP224::new(num)
}
    }

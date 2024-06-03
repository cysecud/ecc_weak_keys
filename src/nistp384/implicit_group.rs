use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};

use super::FieldElement;
use super::scalar384::U384;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitP384{
    pub num:U384,
}
impl Add for ImplicitP384{
    type Output = Self;
    fn add(self:Self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        ImplicitP384{num:sum.mod_fast(Self::prime())}}

}  
impl Mul for ImplicitP384{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        ImplicitP384{num:mul.mod_fast(Self::prime())}}
}
impl Neg for ImplicitP384{
    type Output = Self;
    fn neg(self) -> Self::Output {
        ImplicitP384 {num:Self::prime()-self.num}
    }
}
impl Sub for ImplicitP384{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        
            ImplicitP384 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
        }  
}
impl <'a>FieldElement<'a,U384> for ImplicitP384{
    const PRIME:&'a str ="39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643";
    const PRIME_ROOT:&'a str = "3";
    fn zero()->Self {
        Self { num: U384::zero() }
    }
    fn one()->Self {
        Self { num: U384::one() }
    }
    fn prime()->U384 {
        U384::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U384)->ImplicitP384{
          let prime = U384::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            ImplicitP384 {num:num}
        }
    }
    fn double(&mut self)->ImplicitP384 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U384)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U384::from(2),Self::prime());
        ImplicitP384::new(sqr)
    }
    fn power(&mut self,n: U384) -> Self {
        /*Montgomery ladder */
        
        let bin=U384::to_binary(n);
        let mut p0=ImplicitP384::new(U384::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U384::zero() {
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
    fn inverse (&mut self)-> ImplicitP384 {
        ImplicitP384::power(self, Self::prime()-U384::from(2))
}
fn rand_mod()-> ImplicitP384{
    let prime=U384::from_dec_str(Self::PRIME).expect("error");
    let num=U384::random().mod_fast(prime);
    ImplicitP384::new(num)
}
    }

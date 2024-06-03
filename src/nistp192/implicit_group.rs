use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};
use super::FieldElement;
use super::scalar192::U192;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitP192{
    pub num:U192,
}
impl ImplicitP192{
    pub const DIVISOR_32:[&str; 1]=["191120"];
    //DIVISOR_64 is the same as 32
    pub const DIVISOR_128:[&str;2]=["1828002083835216934513476259280", "656279166350909980926771898430320"];
    //DIVISOR_160 is the same as 128
}

impl Add for ImplicitP192{
    type Output = Self;
    fn add(self:Self,rhs:Self) -> Self::Output{
        let sum=self.num.mod_fast(Self::prime())+rhs.num.mod_fast(Self::prime());
        ImplicitP192{num:sum.mod_fast(Self::prime())}}

}  
impl Mul for ImplicitP192{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mul=self.num.mod_fast(Self::prime())*rhs.num.mod_fast(Self::prime());
        ImplicitP192{num:mul.mod_fast(Self::prime())}}
}
impl Neg for ImplicitP192{
    type Output = Self;
    fn neg(self) -> Self::Output {
        ImplicitP192 {num:Self::prime()-self.num}
    }
}
impl Sub for ImplicitP192{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        
            ImplicitP192 {num: (self.num + other.neg().num).mod_fast(Self::prime())}
        }  
}
impl <'a>FieldElement<'a,U192> for ImplicitP192{
    const PRIME:&'a str ="6277101735386680763835789423176059013767194773182842284081";
    const PRIME_ROOT:&'a str = "3";
    fn zero()->Self {
        Self { num: U192::zero() }
    }
    fn one()->Self {
        Self { num: U192::one() }
    }
    fn prime()->U192 {
        U192::from_dec_str(Self::PRIME).expect("error in Prime Field 192")
    }
    fn new(num:U192)->ImplicitP192{
          let prime = U192::from_dec_str(Self::PRIME).expect("error");
        if num > prime{
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
            } else {
            ImplicitP192 {num:num}
        }
    }
    fn double(&mut self)->ImplicitP192 {
        let doub=self.num.mod_fast(Self::prime())+self.num.mod_fast(Self::prime());
    Self::new(doub.mod_fast(Self::prime()))
    }
    fn multiple(&self,n:U192)->Self {
        let mul=(self.num*n.mod_fast(Self::prime())).mod_fast(Self::prime());
        Self{num:mul}
    }
    fn square(&mut self)->Self{
        let sqr=self.num.pow_mod(U192::from(2),Self::prime());
        ImplicitP192::new(sqr)
    }
    fn power(&mut self,n: U192) -> Self {
        /*Montgomery ladder */
        
        let bin=U192::to_binary(n);
        let mut p0=ImplicitP192::new(U192::one());
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U192::zero() {
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
    fn inverse (&mut self)-> ImplicitP192 {
        ImplicitP192::power(self, Self::prime()-U192::from(2))
}
fn rand_mod()-> ImplicitP192{
    let prime=U192::from_dec_str(Self::PRIME).expect("error");
    let num=U192::random().mod_fast(prime);
    ImplicitP192::new(num)
}
    }

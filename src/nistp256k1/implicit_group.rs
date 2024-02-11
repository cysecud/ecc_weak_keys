use core::ops::{Add, Mul, Sub,Neg};
use core::cmp::{PartialEq,PartialOrd};
use super::FieldElement;
use super::scalar256::U256;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]

pub struct ImplicitP256k1{
    pub num:U256,
}
impl ImplicitP256k1{
    pub const DIVISOR_32:[&str; 1]=["18051648"];
    pub const DIVISOR_64:[&str; 4]=["18051648", "6871154804262114368","10306732206393171552", "15996907278672735013"];
    pub const DIVISOR_128:[&str;8]=["1938057310625759192294976", "3154049060501396119538597952", "9782462315356940942943205974373470912", "41427743093894159295282972951876913728", "192897928780944679218661342807176879546", "225103679124018650938499972290963477356", "257197238374592905624881790409569172728", "300138238832024867917999963054617969808"];
    pub const DIVISOR_160:[&str;2]=["6172733720990229734997162969829660145472", "338624364920977752681389262317185522840540224"];

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

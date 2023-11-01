use crate::affine::AffinePoint;
use crate::projective::ProjectivePoint;
use crate::field::FieldPoint;
pub use crate::ellinit::{Ellinit,CurveParms};
use crate::scalar::U256;

/*Generator of elliptic curve P256r1 */  
pub struct  P256r1{
    pub q:U256,
    pub a:FieldPoint,
    pub b:FieldPoint,
}
impl P256r1 {

pub const XG:&str="55066263022277343669578718895168534326250603453777594175500187360389116729240";
pub const YG:&str="32670510020758816978083085130507043184471273380659243275938904335757337482424";
pub const P:&str="115792089237316195423570985008687907852837564279074904382605163141518161494337";
pub const PRIME_ROOT:&str="7";
    
    pub fn initialize()->Self{ 
        let q:U256=U256::from_dec_str(P256r1::Q).expect("Error in Q");
        let a:FieldPoint=FieldPoint::new(U256::from_dec_str(P256r1::A).expect("Error in coefficient A"),q);
        let b:FieldPoint=FieldPoint::new(U256::from_dec_str(P256r1::B).expect("Error in coefficient B"),q);
        
        Self{a:a,b:b,q:q}
    }
   
}
impl <'a>CurveParms<'a> for P256r1 {
    const Q:&'a str= "115792089210356248762697446949407573530086143415290314195533631308867097853951";
    const A:&'a str= "115792089210356248762697446949407573530086143415290314195533631308867097853948";  //-3 mod q
    const B:&'a str= "41058363725152142129326129780047268409114441015993725554835256314039467401291";

    fn identity(&self)->AffinePoint {
        AffinePoint{ 
            x: FieldPoint{num:U256::zero(),prime:U256::from_dec_str(Self::Q).expect("error in prime Q")}, 
            y: FieldPoint{num:U256::zero(),prime:U256::from_dec_str(Self::Q).expect("error in prime Q")}, 
            infinity: 1 }
    }
    fn generator(&self)->ProjectivePoint 
        {ProjectivePoint {
            x: FieldPoint::new(U256::from_dec_str(Self::XG).expect("Error in x generator"),U256::from_dec_str(Self::Q).expect("error in prime Q x generetor")),
            y: FieldPoint::new(U256::from_dec_str(Self::YG).expect("Error in y generator"),U256::from_dec_str(Self::Q).expect("error in prime Q y generetor")),
            z:FieldPoint::new(U256::one(),U256::from_dec_str(Self::Q).expect("error in prime Q x generetor")),
            infinity: 0,}
    }
    fn zn_prim_root_gen_order(&self)-> U256{
            U256::from_dec_str(Self::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
    fn gen_order(&self)-> U256{
        U256::from_dec_str(Self::P).expect("error in generator's order P!")
    }
    
}

/*
Order of the generator is calculted in pari-gp: ellorder(E,P)
PRIME_ROOT is calculated in pari gp: znprimroot(p).
This is the factorisation of p-1 and can be use to test your key 
with the implicit baby step giant step

The factorisation of p-1 is 

[                           2 4]

[                           3 1]

[                          71 1]

[                         131 1]

[                         373 1]

[                        3407 1]

[                       17449 1]

[                       38189 1]

[                   187019741 1]

[                   622491383 1]

[               1002328039319 1]

[2624747550333869278416773953 1]
*/

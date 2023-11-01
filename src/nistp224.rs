use crate::affine::AffinePoint;
use crate::projective::ProjectivePoint;
use crate::field::FieldPoint;
pub use crate::ellinit::{Ellinit,CurveParms};
use crate::scalar::U256;

/*Structur of elliptic curve p224 */  
pub struct  P224{
    pub q:U256,
    pub a:FieldPoint,
    pub b:FieldPoint,
}
impl P224 {

pub const XG:&str="19277929113566293071110308034699488026831934219452440156649784352033";
pub const YG:&str="19926808758034470970197974370888749184205991990603949537637343198772";
pub const P:&str="26959946667150639794667015087019625940457807714424391721682722368061";
pub const PRIME_ROOT:&str="2";
    
    pub fn initialize()->Self{ 
        let q:U256=U256::from_dec_str(P224::Q).expect("Error in Q");
        let a:FieldPoint=FieldPoint::new(U256::from_dec_str(P224::A).expect("Error in coefficient A"),q);
        let b:FieldPoint=FieldPoint::new(U256::from_dec_str(P224::B).expect("Error in coefficient B"),q);
        
        Self{a:a,b:b,q:q}
    }
   
}
impl <'a>CurveParms<'a> for P224 {
    const Q:&'a str= "26959946667150639794667015087019630673557916260026308143510066298881";
    const A:&'a str="26959946667150639794667015087019630673557916260026308143510066298878";  
    const B:&'a str="18958286285566608000408668544493926415504680968679321075787234672564";

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

[                                                          2 2]

[                                                          3 6]

[                                                          5 1]

[                                                         17 1]

[                                                       2153 1]

[50520606258875818707470860153287666700917696099933389351507 1]


*/

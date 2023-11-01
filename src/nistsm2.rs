use crate::affine::AffinePoint;
use crate::projective::ProjectivePoint;
use crate::field::FieldPoint;
pub use crate::ellinit::{Ellinit,CurveParms};
use crate::scalar::U256;

/*Generator of elliptic curve Sm2 */  
pub struct  Sm2{
    pub q:U256,
    pub a:FieldPoint,
    pub b:FieldPoint,
}
impl Sm2 {

pub const XG:&str="22963146547237050559479531362550074578802567295341616970375194840604139615431";
pub const YG:&str="85132369209828568825618990617112496413088388631904505083283536607588877201568";
pub const P:&str="115792089210356248756420345214020892766061623724957744567843809356293439045923";
pub const PRIME_ROOT:&str="3";
    
    pub fn initialize()->Self{ 
        let q:U256=U256::from_dec_str(Sm2::Q).expect("Error in Q");
        let a:FieldPoint=FieldPoint::new(U256::from_dec_str(Sm2::A).expect("Error in coefficient A"),q);
        let b:FieldPoint=FieldPoint::new(U256::from_dec_str(Sm2::B).expect("Error in coefficient B"),q);
        
        Self{a:a,b:b,q:q}
    }
   
}
impl <'a>CurveParms<'a> for Sm2 {
    const Q:&'a str= "115792089210356248756420345214020892766250353991924191454421193933289684991999";
    const A:&'a str= "115792089210356248756420345214020892766250353991924191454421193933289684991996";  //-3 mod q
    const B:&'a str= "18505919022281880113072981827955639221458448578012075254857346196103069175443";

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

[                                                           2 1]

[                                                           3 1]

[                                                        7759 1]

[                                                       14057 1]

[                                                  1413296869 1]

[125197554539772723432468576818475380947471091418362477724521 1]

*/

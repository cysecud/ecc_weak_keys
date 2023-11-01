
use crate::affine::AffinePoint;
use crate::projective::ProjectivePoint;
use crate::field::FieldPoint;
use crate::ellinit::CurveParms;
use crate::scalar::U256;
/*Structur of elliptic curve p192 */  
pub struct  P192{
    pub q:U256,
    pub a:FieldPoint,
    pub b:FieldPoint,
}
impl P192 {

pub const XG:&str="602046282375688656758213480587526111916698976636884684818";
pub const YG:&str="174050332293622031404857552280219410364023488927386650641";
pub const P:&str="6277101735386680763835789423176059013767194773182842284081";
pub const PRIME_ROOT:&str="3";

pub fn initialize()->Self{ 
    let q:U256=U256::from_dec_str(P192::Q).expect("Error in Q");
    let a:FieldPoint=FieldPoint::new(U256::from_dec_str(P192::A).expect("Error in coefficient A"),q);
    let b:FieldPoint=FieldPoint::new(U256::from_dec_str(P192::B).expect("Error in coefficient B"),q);
    
    Self{a:a,b:b,q:q}
}
} 
impl <'a>CurveParms<'a> for P192 {
    const Q:&'a str="6277101735386680763835789423207666416083908700390324961279"; 
    const A:&'a str = "2455155546008943817740293915197451784769108058161191238065";
    const B:&'a str = "2455155546008943817740293915197451784769108058161191238065";
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

[                           5 1]

[                        2389 1]

[   9564682313913860059195669 1]

[3433859179316188682119986911 1]


*/
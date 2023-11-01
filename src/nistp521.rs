use crate::affine::AffinePoint;
use crate::projective::ProjectivePoint;
use crate::field::FieldPoint;
pub use crate::ellinit::{Ellinit,CurveParms};
use crate::scalar::U256;

/*Generator of elliptic curve P521 */  
pub struct  P521{
    pub q:U256,
    pub a:FieldPoint,
    pub b:FieldPoint,
}
impl P521 {

pub const XG:&str="2661740802050217063228768716723360960729859168756973147706671368418802944996427808491545080627771902352094241225065558662157113545570916814161637315895999846";
pub const YG:&str="3757180025770020463545507224491183603594455134769762486694567779615544477440556316691234405012945539562144444537289428522585666729196580810124344277578376784";
pub const P:&str="6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449";
pub const PRIME_ROOT:&str="3";
    
    pub fn initialize()->Self{ 
        let q:U256=U256::from_dec_str(P521::Q).expect("Error in Q");
        let a:FieldPoint=FieldPoint::new(U256::from_dec_str(P521::A).expect("Error in coefficient A"),q);
        let b:FieldPoint=FieldPoint::new(U256::from_dec_str(P521::B).expect("Error in coefficient B"),q);
        
        Self{a:a,b:b,q:q}
    }
   
}
impl <'a>CurveParms<'a> for P521 {
    const Q:&'a str= "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151";
    const A:&'a str= "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057148";  //-3 mod q
    const B:&'a str= "1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984";

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

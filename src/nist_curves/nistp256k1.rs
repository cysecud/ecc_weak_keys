use crate::affine::AffinePoint;
use crate::projective::ProjectivePoint;
use crate::field::FieldPoint;
use crate::ellinit::Ellinit;
use crate::scalar;

pub type Scalar=scalar::U256;
/*Here we declair ol parameters of p224 elliptic curve */
/* 
pub const A:&str="0";  
pub const B:&str="7";
pub const Q:&str="115792089237316195423570985008687907853269984665640564039457584007908834671663";

        pub const XG:&str="55066263022277343669578718895168534326250603453777594175500187360389116729240";
        pub const YG:&str="32670510020758816978083085130507043184471273380659243275938904335757337482424";
        pub const P:&str="115792089237316195423570985008687907852837564279074904382605163141518161494337";
        pub const PRIME_ROOT:&str="7";
 */
/*Generator of elliptic curve p224 */  
pub struct P256k1;
impl P256k1 {
pub const A:&str="0";  
pub const B:&str="7";
pub const Q:&str="115792089237316195423570985008687907853269984665640564039457584007908834671663";

        pub const XG:&str="55066263022277343669578718895168534326250603453777594175500187360389116729240";
        pub const YG:&str="32670510020758816978083085130507043184471273380659243275938904335757337482424";
        pub const P:&str="115792089237316195423570985008687907852837564279074904382605163141518161494337";
        pub const PRIME_ROOT:&str="7";
    
    pub fn initialize(self)->Ellinit{ 
        let q:Scalar=Scalar::from_dec_str(P256k1::Q).expect("Error in Q");
        let a:FieldPoint=FieldPoint::new(Scalar::from_dec_str(P256k1::A).expect("Error in coefficient A"),q);
        let b:FieldPoint=FieldPoint::new(Scalar::from_dec_str(P256k1::B).expect("Error in coefficient B"),q);
        
        Ellinit::new(a,b,q)
    }
}
impl Ellinit for P256k1 {
     /* pub fn identity(self)->AffinePoint {
        AffinePoint{ 
            x: FieldPoint{num:Scalar::zero(),prime:Scalar::from_dec_str(P256k1::Q).expect("error in prime Q")}, 
            y: FieldPoint{num:Scalar::zero(),prime:Scalar::from_dec_str(P256k1::Q).expect("error in prime Q")}, 
            infinity: 1 }
    } */
    fn generator(self)->ProjectivePoint 
        {ProjectivePoint {
            x: FieldPoint::new(Scalar::from_dec_str(P256k1::XG).expect("Error in x generator"),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q x generetor")),
            y: FieldPoint::new(Scalar::from_dec_str(P256k1::YG).expect("Error in y generator"),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q y generetor")),
            z:FieldPoint::new(Scalar::one(),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q x generetor")),
            infinity: 0,}
    }
     fn zn_prim_root_gen_order(self)-> Scalar{
            Scalar::from_dec_str(P256k1::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
     fn gen_order(self)-> Scalar{
        Scalar::from_dec_str(P256k1::P).expect("error in generator's order P!")
    }
    
}
impl  ProjectivePoint {
    pub fn identity()->ProjectivePoint {
       Self {
            x: FieldPoint::new(Scalar::zero(),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q x projective identity")),
            y: FieldPoint::new(Scalar::zero(),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q z projective identiy")),
            z: FieldPoint::new(Scalar::zero(),Scalar::from_dec_str(P256k1::Q).expect("error in prime Qz projective identity")),
            infinity: 1,}
    }
}
impl  AffinePoint {
    pub fn identity()->AffinePoint {
       Self {
            x: FieldPoint::new(Scalar::zero(),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q x projective identity")),
            y: FieldPoint::new(Scalar::zero(),Scalar::from_dec_str(P256k1::Q).expect("error in prime Q z projective identiy")),
            infinity: 1,}
    }
}

/*Order of the generator. It can be calculted in pari-gp: ellorder(E,P) */
/*PRIME_ROOT can be calculated in pari gp. It is a primitive root of the base field.
and we need it for the implicit baby step giant step

The factorisation of p-1 is 

[                                                          2 2]

[                                                          3 6]

[                                                          5 1]

[                                                         17 1]

[                                                       2153 1]

[50520606258875818707470860153287666700917696099933389351507 1]

*/

pub mod scalar256;
pub mod field_p256k1;
pub mod implicit_group;

use self::scalar256::{U256,MathResult};
use self::field_p256k1::FieldP256k1;
use self::implicit_group::ImplicitP256k1;
use std::ops::Neg;

use crate::affine::AffinePoint;
use crate::FieldElement;
use crate::ellinit::{CurveParms, KeyPairs};
use crate::ellinit::Curve;
use crate::projective::ProjectivePoint;


#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
/*Structur of elliptic curve p192 */  
pub struct  P256k1{
    pub q:U256,
    pub a:FieldP256k1,
    pub b:FieldP256k1,
}

impl <'a>CurveParms<'a> for P256k1 {
    const A:&'a str = "0";
    const B:&'a str = "7";
    const XG:&'a str= "55066263022277343669578718895168534326250603453777594175500187360389116729240";
    const YG:&'a str= "32670510020758816978083085130507043184471273380659243275938904335757337482424";

}
impl Curve<FieldP256k1,U256,MathResult> for P256k1 {
    fn initialize()->Self {
        let q=U256::from_dec_str(FieldP256k1::PRIME).expect("error in Prime FieldP256k1 inizializing U256 curve");
        let a=FieldP256k1::new(U256::from_dec_str(Self::A).expect("Error in coefficient A"));
        let b=FieldP256k1::new(U256::from_dec_str(Self::B).expect("Error in coefficient B"));
            
        Self{a,b,q}
    }
    fn identity(&self)->AffinePoint<FieldP256k1>{
        AffinePoint{ 
            x: FieldP256k1::zero(), 
            y: FieldP256k1::one(), 
            infinity: 1 }
    }
    fn generator(&self)->ProjectivePoint<FieldP256k1> 
        {ProjectivePoint {
            x: FieldP256k1::new(U256::from_dec_str(Self::XG).expect("Error in x generator")),
            y: FieldP256k1::new(U256::from_dec_str(Self::YG).expect("Error in y generator")),
            z:FieldP256k1::one(),
            infinity: 0,}
    }
    fn ellnegate(&self,p:AffinePoint<FieldP256k1>)->AffinePoint<FieldP256k1> {
        AffinePoint { x: p.x, y: p.y.neg(), infinity: p.infinity }
    }
    fn ellnegate_proj(&self,p:ProjectivePoint<FieldP256k1>)->ProjectivePoint<FieldP256k1> {
        ProjectivePoint { x: p.x, y: p.y.neg(), z:p.z,infinity: p.infinity }
    }
     fn zn_prim_root_gen_order(&self)-> U256{
            U256::from_dec_str(ImplicitP256k1::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
     fn gen_order(&self)-> U256{
        U256::from_dec_str(ImplicitP256k1::PRIME).expect("error in generator's order P!")
    }

    fn ellisoncurve(&self,mut p:AffinePoint<FieldP256k1>)->bool {
        let a=FieldP256k1::new(U256::from_dec_str(Self::A).expect("error in A"));
        let b=FieldP256k1::new(U256::from_dec_str(Self::B).expect("error in B"));
        
        let x3=p.x.power(U256::from(3));
        let ax= a*p.x;
        let y2=x3+ax+b;
        let check:Option<U256>=U256::check_sqrt_mod_prime(y2.num, FieldP256k1::prime());
        if check.is_some() {
            if p.y.num==check.unwrap() || p.y.num==(FieldP256k1::prime()-check.unwrap()) {true} else{false}
        }else{false}
    }

    fn random (&self)->AffinePoint<FieldP256k1> {
        let mut x:FieldP256k1;
        let y:FieldP256k1;
        let a:FieldP256k1=FieldP256k1::new(U256::from_dec_str(Self::A).expect("error in A"));
        let b:FieldP256k1=FieldP256k1::new(U256::from_dec_str(Self::B).expect("error in B"));
        x=loop{
            x=FieldP256k1::rand_mod();
            let x3=FieldP256k1::power(&mut x,U256::from(3));
            let ax=a*x;
            let y2=x3+ax+b;
            let r:Option<U256>=U256::check_sqrt_mod_prime(y2.num,FieldP256k1::prime());
            if r.is_some() 
                {y=FieldP256k1::new(r.unwrap());
                    break x;}
        };
        AffinePoint::new(x,y)    }

    fn private_key(&self)->U256 {
        FieldP256k1::rand_mod().num    }

    fn key_pairs(&self)->KeyPairs<U256,AffinePoint<FieldP256k1>> {
            let sk=self.private_key();
            let pk = self.ellmul(&mut self.to_affine(self.generator()), sk);
            KeyPairs{sk,pk}
    }
    fn to_affine(&self,mut p:ProjectivePoint<FieldP256k1>)->AffinePoint<FieldP256k1> {
        let x_aff=p.x*FieldP256k1::inverse( &mut p.z);
        let y_aff=p.y*FieldP256k1::inverse(&mut p.z);
        AffinePoint::new(x_aff,y_aff)
    }

    fn jc_to_affine(&self,mut p:ProjectivePoint<FieldP256k1>)->AffinePoint<FieldP256k1> {
        let mut inv_z=p.z.inverse();
        let z=inv_z.square();
        let t=z*inv_z;
        let aff_x=p.x*z;
        let aff_y=p.y*t;
    AffinePoint::new(aff_x,aff_y)
    }

    fn proj_identity(&self)->ProjectivePoint<FieldP256k1> {
        ProjectivePoint{ 
            x: FieldP256k1::zero(), 
            y: FieldP256k1::one(),
            z: FieldP256k1::zero(), 
            infinity: 1 }
    }
    fn proj_elladd_jc(&self,mut p:ProjectivePoint<FieldP256k1>,mut q:ProjectivePoint<FieldP256k1>)->ProjectivePoint<FieldP256k1> {
            /*Jacobian and Chudnovsky Jacobian coordinates
        With Jacobian coordinates the curve E is given by
                    Y^2 = X^3 + a4XZ 4 + a6Z^6 .
        The point (X1 : Y1 : Z1 ) on E corresponds to the affine point (X1 /Z1^2 , Y1 /Z1^3 ) when Z1 is not 0 and to
        the point at infinity otherwise. The opposite of (X1 : Y1 : Z1 ) is (X1 : −Y1 : Z1 ). 
        The complexities are 12M + 4S for an addition and 4M + 6S for a doubling. If one of the points is
        given in the form (X1 : Y1 : 1) the costs for addition reduce to 8M + 3S.*/
        if p==self.proj_identity() {return q;} 
        else if q==self.proj_identity() {return p;}
        else if p==self.ellnegate_proj(q) {return self.proj_identity();}
        else{
            let a=p.x*q.z.square();
            let b=q.x*p.z.square();
            let c=p.y*q.z.square()*q.z;
            let d=q.y*p.z.square()*p.z;
            let mut e=b-a;
            let mut f=d-c;
            let x3=-(e.square()*e)-(a*e.square().double())+f.square();
            let y3=-(c*e.square()*e)+f*(a*e.square()-x3);
            let z3=p.z*q.z*e;
            ProjectivePoint::new(x3,y3,z3)
        }
    }
    fn elladd(&self,aff_p:AffinePoint<FieldP256k1>,aff_q: AffinePoint<FieldP256k1>)->AffinePoint<FieldP256k1> {
            //let a4=FieldP256k1::new(U256::from_dec_str(Self::A).expect("error in A elladd"),aff_p.x.prime);
            //let q=aff_p.x.prime;
            if aff_p==Self::identity(&self) {return aff_q;}
            else if aff_q==Self::identity(self) {return aff_p;}
            else if aff_p==self.ellnegate(aff_q) {return  self.identity();}
    /*         else if q-U256::from(3)==a4.num {self.to_affine(self.proj_elladd_3(proj_p, proj_q))}
     */        else{
            let proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP256k1::one());
            let proj_q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldP256k1::one());
            self.jc_to_affine(self.proj_elladd_jc(proj_p, proj_q))}
    
}
fn proj_elldouble_jc(&self, proj_p:&mut ProjectivePoint<FieldP256k1>)->ProjectivePoint<FieldP256k1> {
    
    let a4:FieldP256k1=FieldP256k1::new(U256::from_dec_str(Self::A).expect("error in A proj_elldouble_jc"));

    let mut a=FieldP256k1::multiple(&(proj_p.x*(proj_p.y.square())),U256::from(4));
    let mut b:FieldP256k1;
    if FieldP256k1::prime()-U256::from(3)==a4.num {b=FieldP256k1::multiple(&((proj_p.x-proj_p.z.square())*(proj_p.x+proj_p.z.square())),U256::from(3));}
    else if a4.num==U256::zero(){b=FieldP256k1::multiple(&(proj_p.x.square()), U256::from(3))}
    else {b=FieldP256k1::multiple(&proj_p.x.square(),U256::from(3))+a4*proj_p.z.power(U256::from(4));}
    let x3=b.square()-a.double();
    let y3=-(proj_p.y.power(U256::from(4))).multiple(U256::from(8))+b*(a-x3);
    let z3=proj_p.y*proj_p.z.double();
    ProjectivePoint ::new(x3,y3, z3)
}
fn elldouble(&self,aff_p:&mut AffinePoint<FieldP256k1>)->AffinePoint<FieldP256k1> {
    let mut proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP256k1::one());
    self.jc_to_affine(self.proj_elldouble_jc(&mut proj_p))
}
fn ellmul(&self,p:&mut AffinePoint<FieldP256k1>,k:U256)->AffinePoint<FieldP256k1> {
    {

        /**************Montgomery ladder******************
[https://www.matthieurivain.com/files/jcen11b.pdf] 
Algorithm 3. */

let mut p0=self.proj_identity();
let mut p1= ProjectivePoint::new(p.x,p.y,FieldP256k1::one());
let bin=U256::to_binary(k);
for i in bin.iter(){
    if *i==U256::zero(){
        p1=self.proj_elladd_jc(p0, p1);

        p0=self.proj_elldouble_jc( &mut p0);}
    else {
        p0=self.proj_elladd_jc(p0, p1);
        p1=Self::proj_elldouble_jc(&self,&mut p1);}
    }
    self.jc_to_affine(p0)
}   

}
fn bsgs(&self,q:&mut AffinePoint<FieldP256k1>, d:U256)->MathResult {
        //This is the order of the generator point gen, it is calculated in pari-gp
	let ord:U256 = self.gen_order();
	//let mut z = FieldP256k1::znprimroot(&ord); /*generatore di Fp* dove p=ord(E,P)*/
	//This is a primitive root mod p, it is calculated in pari-gp
	let mut z= ImplicitP256k1::new(U256::from(self.zn_prim_root_gen_order()));
	let mut zd = ImplicitP256k1::power(&mut z,(ord-U256::one())/d);
	//calculate the inverse of the d-order generator
    let mut inv_zd=zd.inverse();
	let mut res:MathResult=None;
    let mut bs: Vec<(u128,AffinePoint<FieldP256k1>)>=Vec::new();
    let mut gs: Vec<(u128,AffinePoint<FieldP256k1>)>=Vec::new();
	let  m =U256::integer_sqrt(&d)+U256::one();//ceil(d)
	'outer: for i in 0..m.as_u128() {
        bs.push((i,self.ellmul(&mut self.jc_to_affine(self.generator()), ImplicitP256k1::power(&mut zd,U256::from(i)).num)));
        bs.sort_by_key(|key: &(u128, AffinePoint<FieldP256k1>)| key.1 );
		gs.push((i,self.ellmul(q,ImplicitP256k1::power(&mut inv_zd,m*U256::from(i)).num)));
        for j in 0..bs.len(){
			let search=bs.binary_search_by_key(&gs[j].1,|&(_a,b)|b);
            if search.is_ok()
			{
            let i =gs[j].0;
            let j=bs[search.unwrap()].0;
            println!("Baby Step Giant Step found a match:i={},j={}",i,j);
			let alpha=ImplicitP256k1::power(&mut zd,m*U256::from(i)+U256::from(j)%d);
			println!("The key is weak. Your secret key is {}",alpha.num);
            res=Some(alpha.num);
			break 'outer;
        }}}
			if res.is_some() {return res} else {return None;}	
}          

}
impl P256k1 {
    pub fn test_key(&self, q:&mut AffinePoint<FieldP256k1>,bound:usize)->Option<U256>{
        let div :Vec<U256>=match bound {
            32usize|64usize|128usize|160usize => ImplicitP256k1::DIVISOR_32.iter().map(|x| U256::from_dec_str(x).expect("error in divisor 32")).collect(),
            _=>panic!("bound must be a usize 32, 64,128 or 160!")
            
        };
        println!("div is {:?}",div);
        for item in div {
            let alfa=self.bsgs(q, item);
            if alfa.is_some() {return alfa}
        }
        None
        
    }
}



pub mod scalar256;
pub mod field_sm256;
pub mod implicit_group;

use self::scalar256::{U256,MathResult};
use self::field_sm256::FieldSM256;
use self::implicit_group::ImplicitSM256;
use std::ops::Neg;

use crate::affine::AffinePoint;
use crate::FieldElement;
use crate::ellinit::{CurveParms, KeyPairs};
use crate::ellinit::Curve;
use crate::projective::ProjectivePoint;


#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
/*Structur of elliptic curve p192 */  
pub struct  SM256{
    pub q:U256,
    pub a:FieldSM256,
    pub b:FieldSM256,
}

impl <'a>CurveParms<'a> for SM256 {
    const A:&'a str = "115792089210356248756420345214020892766250353991924191454421193933289684991996";
    const B:&'a str = "18505919022281880113072981827955639221458448578012075254857346196103069175443";
    const XG:&'a str= "22963146547237050559479531362550074578802567295341616970375194840604139615431";
    const YG:&'a str= "85132369209828568825618990617112496413088388631904505083283536607588877201568";

}
impl Curve<FieldSM256,U256,MathResult> for SM256 {
    fn initialize()->Self {
        let q=U256::from_dec_str(FieldSM256::PRIME).expect("error in Prime FieldSM256 inizializing U256 curve");
        let a=FieldSM256::new(U256::from_dec_str(Self::A).expect("Error in coefficient A"));
        let b=FieldSM256::new(U256::from_dec_str(Self::B).expect("Error in coefficient B"));
            
        Self{a:a,b:b,q:q}
    }
    fn identity(&self)->AffinePoint<FieldSM256>{
        AffinePoint{ 
            x: FieldSM256::zero(), 
            y: FieldSM256::one(), 
            infinity: 1 }
    }
    fn generator(&self)->ProjectivePoint<FieldSM256> 
        {ProjectivePoint {
            x: FieldSM256::new(U256::from_dec_str(Self::XG).expect("Error in x generator")),
            y: FieldSM256::new(U256::from_dec_str(Self::YG).expect("Error in y generator")),
            z:FieldSM256::one(),
            infinity: 0,}
    }
    fn ellnegate(&self,p:AffinePoint<FieldSM256>)->AffinePoint<FieldSM256> {
        AffinePoint { x: p.x, y: p.y.neg(), infinity: p.infinity }
    }
    fn ellnegate_proj(&self,p:ProjectivePoint<FieldSM256>)->ProjectivePoint<FieldSM256> {
        ProjectivePoint { x: p.x, y: p.y.neg(), z:p.z,infinity: p.infinity }
    }
     fn zn_prim_root_gen_order(&self)-> U256{
            U256::from_dec_str(ImplicitSM256::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
     fn gen_order(&self)-> U256{
        U256::from_dec_str(ImplicitSM256::PRIME).expect("error in generator's order P!")
    }

    fn ellisoncurve(&self,mut p:AffinePoint<FieldSM256>)->bool {
/*         let q:U256=U256::from_dec_str(Self::Q).expect("error in Q random");
 */     let a=FieldSM256::new(U256::from_dec_str(Self::A).expect("error in A"));
        let b=FieldSM256::new(U256::from_dec_str(Self::B).expect("error in B"));
        
        let x3=p.x.power(U256::from(3));
        let ax= a*p.x;
        let y2=x3+ax+b;
        let check:Option<U256>=U256::check_sqrt_mod_prime(y2.num, FieldSM256::prime());
        if check.is_some() {
            if p.y.num==check.unwrap() || p.y.num==(FieldSM256::prime()-check.unwrap()) {true} else{false}
        }else{false}
    }

    fn random (&self)->AffinePoint<FieldSM256> {
        let mut x:FieldSM256;
        let y:FieldSM256;
        let a:FieldSM256=FieldSM256::new(U256::from_dec_str(Self::A).expect("error in A"));
        let b:FieldSM256=FieldSM256::new(U256::from_dec_str(Self::B).expect("error in B"));
        x=loop{
            x=FieldSM256::rand_mod();
            let x3=FieldSM256::power(&mut x,U256::from(3));
            let ax=a*x;
            let y2=x3+ax+b;
            let r:Option<U256>=U256::check_sqrt_mod_prime(y2.num,FieldSM256::prime());
            if r.is_some() 
                {y=FieldSM256::new(r.unwrap());
                    break x;}
        };
        AffinePoint::new(x,y)    }

    fn private_key(&self)->U256 {
        FieldSM256::rand_mod().num    }

    fn key_pairs(&self)->KeyPairs<U256,AffinePoint<FieldSM256>> {
        let sk=self.private_key();
        let pk = self.ellmul(&mut self.to_affine(self.generator()), sk);
        KeyPairs{sk,pk}
    }

    fn to_affine(&self,mut p:ProjectivePoint<FieldSM256>)->AffinePoint<FieldSM256> {
        let x_aff=p.x*FieldSM256::inverse( &mut p.z);
        let y_aff=p.y*FieldSM256::inverse(&mut p.z);
        AffinePoint::new(x_aff,y_aff)
    }

    fn jc_to_affine(&self,mut p:ProjectivePoint<FieldSM256>)->AffinePoint<FieldSM256> {
        let mut inv_z=p.z.inverse();
        let z=inv_z.square();
        let t=z*inv_z;
        let aff_x=p.x*z;
        let aff_y=p.y*t;
    AffinePoint::new(aff_x,aff_y)
    }

    fn proj_identity(&self)->ProjectivePoint<FieldSM256> {
        ProjectivePoint{ 
            x: FieldSM256::zero(), 
            y: FieldSM256::one(),
            z: FieldSM256::zero(), 
            infinity: 1 }
    }
    fn proj_elladd_jc(&self,mut p:ProjectivePoint<FieldSM256>,mut q:ProjectivePoint<FieldSM256>)->ProjectivePoint<FieldSM256> {
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
    fn elladd(&self,aff_p:AffinePoint<FieldSM256>,aff_q: AffinePoint<FieldSM256>)->AffinePoint<FieldSM256> {
            //let a4=FieldSM256::new(U256::from_dec_str(Self::A).expect("error in A elladd"),aff_p.x.prime);
            //let q=aff_p.x.prime;
            if aff_p==Self::identity(&self) {return aff_q;}
            else if aff_q==Self::identity(self) {return aff_p;}
            else if aff_p==self.ellnegate(aff_q) {return  self.identity();}
    /*         else if q-U256::from(3)==a4.num {self.to_affine(self.proj_elladd_3(proj_p, proj_q))}
     */        else{
            let proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldSM256::one());
            let proj_q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldSM256::one());
            self.jc_to_affine(self.proj_elladd_jc(proj_p, proj_q))}
    
}
fn proj_elldouble_jc(&self, proj_p:&mut ProjectivePoint<FieldSM256>)->ProjectivePoint<FieldSM256> {
    
    let a4:FieldSM256=FieldSM256::new(U256::from_dec_str(Self::A).expect("error in A proj_elldouble_jc"));

    let mut a=FieldSM256::multiple(&(proj_p.x*(proj_p.y.square())),U256::from(4));
    let mut b:FieldSM256;
    if FieldSM256::prime()-U256::from(3)==a4.num {b=FieldSM256::multiple(&((proj_p.x-proj_p.z.square())*(proj_p.x+proj_p.z.square())),U256::from(3));}
    else if a4.num==U256::zero(){b=FieldSM256::multiple(&(proj_p.x.square()), U256::from(3))}
    else {b=FieldSM256::multiple(&proj_p.x.square(),U256::from(3))+a4*proj_p.z.power(U256::from(4));}
    let x3=b.square()-a.double();
    let y3=-(proj_p.y.power(U256::from(4))).multiple(U256::from(8))+b*(a-x3);
    let z3=proj_p.y*proj_p.z.double();
    ProjectivePoint ::new(x3,y3, z3)
}
fn elldouble(&self,aff_p:&mut AffinePoint<FieldSM256>)->AffinePoint<FieldSM256> {
    let mut proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldSM256::one());
    self.jc_to_affine(self.proj_elldouble_jc(&mut proj_p))
}
fn ellmul(&self,p:&mut AffinePoint<FieldSM256>,k:U256)->AffinePoint<FieldSM256> {
    {

        /**************Montgomery ladder******************
[https://www.matthieurivain.com/files/jcen11b.pdf] 
Algorithm 3. */

let mut p0=self.proj_identity();
let mut p1= ProjectivePoint::new(p.x,p.y,FieldSM256::one());
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
fn bsgs(&self, q:&mut AffinePoint<FieldSM256>, d:U256)->MathResult {
        //This is the order of the generator point gen, it is calculated in pari-gp
	let ord:U256 = self.gen_order();
	//let mut z = FieldSM256::znprimroot(&ord); /*generatore di Fp* dove p=ord(E,P)*/
	//This is a primitive root mod p, it is calculated in pari-gp
	let mut z= ImplicitSM256::new(U256::from(self.zn_prim_root_gen_order()));
	let mut zd = ImplicitSM256::power(&mut z,(ord-U256::one())/d);
	//calculate the inverse of the d-order generator
    let mut inv_zd=zd.inverse();
	let mut res:MathResult=None;
    let mut bs: Vec<(u128,AffinePoint<FieldSM256>)>=Vec::new();
    let mut gs: Vec<(u128,AffinePoint<FieldSM256>)>=Vec::new();
	let  m =U256::integer_sqrt(&d)+U256::one();//ceil(d)
	'outer: for i in 0..m.as_u128() {
        bs.push((i,self.ellmul(&mut self.jc_to_affine(self.generator()), ImplicitSM256::power(&mut zd,U256::from(i)).num)));
        bs.sort_by_key(|key: &(u128, AffinePoint<FieldSM256>)| key.1 );
		gs.push((i,self.ellmul(q,ImplicitSM256::power(&mut inv_zd,m*U256::from(i)).num)));
        for j in 0..bs.len(){
			let search=bs.binary_search_by_key(&gs[j].1,|&(_a,b)|b);
            if search.is_ok()
			{
            let i =gs[j].0;
            let j=bs[search.unwrap()].0;
            println!("Baby Step Giant Step found a match:i={},j={}",i,j);
			let alpha=ImplicitSM256::power(&mut zd,m*U256::from(i)+U256::from(j)%d);
			println!("The key is weak. Your secret key is {}",alpha.num);
            res=Some(alpha.num);
			break 'outer;
        }}}
			if res.is_some() {return res} else {return None;}	
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

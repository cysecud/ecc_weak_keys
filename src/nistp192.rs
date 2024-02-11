
pub mod scalar192;
pub mod field_p192;
pub mod implicit_group;

use self::scalar192::{U192,MathResult};
use self::field_p192::FieldP192;
use self::implicit_group::ImplicitP192;
use std::ops::Neg;
use std::str::FromStr;

use crate::FieldElement;
use crate::ellinit::CurveParms;
use crate::ellinit::Curve;
use crate::projective::ProjectivePoint;
use crate::affine::AffinePoint;


/*Structur of elliptic curve p192 */  
pub struct  P192{
    pub q:U192,
    pub a:FieldP192,
    pub b:FieldP192,
}
impl <'a>CurveParms<'a> for P192 {
/*     const Q:&'a str= "6277101735386680763835789423207666416083908700390324961279"; 
 */    const A:&'a str = "6277101735386680763835789423207666416083908700390324961276";
    const B:&'a str = "2455155546008943817740293915197451784769108058161191238065";
    const XG:&'a str= "602046282375688656758213480587526111916698976636884684818";
    const YG:&'a str= "174050332293622031404857552280219410364023488927386650641";
    const P:&'a str= "6277101735386680763835789423176059013767194773182842284081";
    const PRIME_ROOT:&'a str="3";
}
impl Curve<FieldP192,U192,MathResult> for P192 {
    fn initialize()->Self {
        let q=U192::from_dec_str(FieldP192::PRIME).expect("error in Prime FieldP192 inizializing P192 curve");
        let a=FieldP192::new(U192::from_dec_str(Self::A).expect("Error in coefficient A"));
        let b=FieldP192::new(U192::from_dec_str(Self::B).expect("Error in coefficient B"));
            
        Self{a:a,b:b,q:q}
    }
    fn identity(&self)->AffinePoint<FieldP192>{
        AffinePoint{ 
            x: FieldP192::zero(), 
            y: FieldP192::one(), 
            infinity: 1 }
    }
    fn generator(&self)->ProjectivePoint<FieldP192> 
        {ProjectivePoint {
            x: FieldP192::new(U192::from_dec_str(Self::XG).expect("Error in x generator")),
            y: FieldP192::new(U192::from_dec_str(Self::YG).expect("Error in y generator")),
            z:FieldP192::one(),
            infinity: 0,}
    }
    fn ellnegate(&self,p:AffinePoint<FieldP192>)->AffinePoint<FieldP192> {
        AffinePoint { x: p.x, y: p.y.neg(), infinity: p.infinity }
    }
    fn ellnegate_proj(&self,p:ProjectivePoint<FieldP192>)->ProjectivePoint<FieldP192> {
        ProjectivePoint { x: p.x, y: p.y.neg(), z:p.z,infinity: p.infinity }
    }
     fn zn_prim_root_gen_order(&self)-> U192{
            U192::from_dec_str(Self::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
     fn gen_order(&self)-> U192{
        U192::from_dec_str(Self::P).expect("error in generator's order P!")
    }

    fn ellisoncurve(&self,mut p:AffinePoint<FieldP192>)->bool {
/*         let q:U192=U192::from_dec_str(Self::Q).expect("error in Q random");
 */     let a=FieldP192::new(U192::from_dec_str(Self::A).expect("error in A"));
        let b=FieldP192::new(U192::from_dec_str(Self::B).expect("error in B"));
        
        let x3=p.x.power(U192::from(3));
        let ax= a*p.x;
        let y2=x3+ax+b;
        let check:Option<U192>=U192::check_sqrt_mod_prime(y2.num, FieldP192::prime());
        if check.is_some() {
            if p.y.num==check.unwrap() || p.y.num==(FieldP192::prime()-check.unwrap()) {true} else{false}
        }else{false}
    }

    fn random (&self)->AffinePoint<FieldP192> {
        let mut x:FieldP192;
        let y:FieldP192;
        let a:FieldP192=FieldP192::new(U192::from_dec_str(Self::A).expect("error in A"));
        let b:FieldP192=FieldP192::new(U192::from_dec_str(Self::B).expect("error in B"));
        x=loop{
            x=FieldP192::rand_mod();
            let x3=FieldP192::power(&mut x,U192::from(3));
            let ax=a*x;
            let y2=x3+ax+b;
            let r:Option<U192>=U192::check_sqrt_mod_prime(y2.num,FieldP192::prime());
            if r.is_some() 
                {y=FieldP192::new(r.unwrap());
                    break x;}
        };
        AffinePoint::new(x,y)    }

    fn private_key(&self)->U192 {
        FieldP192::rand_mod().num    }

    fn public_key(&self,private_key:U192)->AffinePoint<FieldP192> {
        todo!()
    }

    fn to_affine(&self,mut p:ProjectivePoint<FieldP192>)->AffinePoint<FieldP192> {
        let x_aff=p.x*FieldP192::inverse( &mut p.z);
        let y_aff=p.y*FieldP192::inverse(&mut p.z);
        AffinePoint::new(x_aff,y_aff)
    }

    fn jc_to_affine(&self,mut p:ProjectivePoint<FieldP192>)->AffinePoint<FieldP192> {
        let mut inv_z=p.z.inverse();
        let z=inv_z.square();
        let t=z*inv_z;
        let aff_x=p.x*z;
        let aff_y=p.y*t;
    AffinePoint::new(aff_x,aff_y)
    }

    fn proj_identity(&self)->ProjectivePoint<FieldP192> {
        ProjectivePoint{ 
            x: FieldP192::zero(), 
            y: FieldP192::one(),
            z: FieldP192::zero(), 
            infinity: 1 }
    }
    fn proj_elladd_jc(&self,mut p:ProjectivePoint<FieldP192>,mut q:ProjectivePoint<FieldP192>)->ProjectivePoint<FieldP192> {
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
    fn elladd(&self,aff_p:AffinePoint<FieldP192>,aff_q: AffinePoint<FieldP192>)->AffinePoint<FieldP192> {
            //let a4=FieldP192::new(U192::from_dec_str(Self::A).expect("error in A elladd"),aff_p.x.prime);
            //let q=aff_p.x.prime;
            if aff_p==Self::identity(&self) {return aff_q;}
            else if aff_q==Self::identity(self) {return aff_p;}
            else if aff_p==self.ellnegate(aff_q) {return  self.identity();}
    /*         else if q-U192::from(3)==a4.num {self.to_affine(self.proj_elladd_3(proj_p, proj_q))}
     */        else{
            let proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP192::one());
            let proj_q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldP192::one());
            self.jc_to_affine(self.proj_elladd_jc(proj_p, proj_q))}
    
}
fn proj_elldouble_jc(&self, proj_p:&mut ProjectivePoint<FieldP192>)->ProjectivePoint<FieldP192> {
    
    let a4:FieldP192=FieldP192::new(U192::from_dec_str(Self::A).expect("error in A proj_elldouble_jc"));

    let mut a=FieldP192::multiple(&(proj_p.x*(proj_p.y.square())),U192::from(4));
    let mut b:FieldP192;
    if FieldP192::prime()-U192::from(3)==a4.num {b=FieldP192::multiple(&((proj_p.x-proj_p.z.square())*(proj_p.x+proj_p.z.square())),U192::from(3));}
    else if a4.num==U192::zero(){b=FieldP192::multiple(&(proj_p.x.square()), U192::from(3))}
    else {b=FieldP192::multiple(&proj_p.x.square(),U192::from(3))+a4*proj_p.z.power(U192::from(4));}
    let x3=b.square()-a.double();
    let y3=-(proj_p.y.power(U192::from(4))).multiple(U192::from(8))+b*(a-x3);
    let z3=proj_p.y*proj_p.z.double();
    ProjectivePoint ::new(x3,y3, z3)
}
fn elldouble(&self,aff_p:&mut AffinePoint<FieldP192>)->AffinePoint<FieldP192> {
    let mut proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP192::one());
    self.jc_to_affine(self.proj_elldouble_jc(&mut proj_p))
}
fn ellmul(&self,p:&mut AffinePoint<FieldP192>,k:U192)->AffinePoint<FieldP192> {
    {

        /**************Montgomery ladder******************
[https://www.matthieurivain.com/files/jcen11b.pdf] 
Algorithm 3. */

let mut p0=self.proj_identity();
let mut p1= ProjectivePoint::new(p.x,p.y,FieldP192::one());
let bin=U192::to_binary(k);
for i in bin.iter(){
    if *i==U192::zero(){
        p1=self.proj_elladd_jc(p0, p1);

        p0=self.proj_elldouble_jc( &mut p0);}
    else {
        p0=self.proj_elladd_jc(p0, p1);
        p1=Self::proj_elldouble_jc(&self,&mut p1);}
    }
    self.jc_to_affine(p0)
}   

}
fn bsgs(&self,q:&mut AffinePoint<FieldP192>, d:U192)->Option<U192> {
        //This is the order of the generator point gen, it is calculated in pari-gp
	let ord:U192 = self.gen_order();
	//let mut z = FieldP192::znprimroot(&ord); /*generatore di Fp* dove p=ord(E,P)*/
	//This is a primitive root mod p, it is calculated in pari-gp
	let mut z= ImplicitP192::new(U192::from(self.zn_prim_root_gen_order()));
	let mut zd = ImplicitP192::power(&mut z,(ord-U192::one())/d);
	//calculate the inverse of the d-order generator
    let mut inv_zd=zd.inverse();
	let mut res:MathResult=None;
    let mut bs: Vec<(u128,AffinePoint<FieldP192>)>=Vec::new();
    let mut gs: Vec<(u128,AffinePoint<FieldP192>)>=Vec::new();
	let  m =U192::integer_sqrt(&d)+U192::one();//ceil(d)
	'outer: for i in 0..m.as_u128() {
        bs.push((i,self.ellmul(&mut self.jc_to_affine(self.generator()), ImplicitP192::power(&mut zd,U192::from(i)).num)));
        bs.sort_by_key(|key: &(u128, AffinePoint<FieldP192>)| key.1 );
		gs.push((i,self.ellmul(q,ImplicitP192::power(&mut inv_zd,m*U192::from(i)).num)));
        for j in 0..bs.len(){
			let search=bs.binary_search_by_key(&gs[j].1,|&(_a,b)|b);
            if search.is_ok()
			{
               println!("bs is {:?}",bs[search.unwrap()]);
                println!("gs is {:?}",gs[j]);
            let i =gs[j].0;
            let j=bs[search.unwrap()].0;
            println!("Baby Step Giant Step found a match:i={},j={}",i,j);
			let alpha=ImplicitP192::power(&mut zd,m*U192::from(i)+U192::from(j)%d);
			println!("The key is weak. Your secret key is {}",alpha.num);
            res=Some(alpha.num);
			break 'outer;
        }}}
			if res.is_some() {return res} else {return None;}	
}          

}
impl P192 {
    pub fn test_key(&self, q:&mut AffinePoint<FieldP192>,bound:usize)->Option<U192>{
        let div :Vec<U192>=match bound {
            32usize => ImplicitP192::DIVISOR_32.iter().map(|x| U192::from_dec_str(x).expect("error in divisor 32")).collect(),
            128usize=>ImplicitP192::DIVISOR_128.iter().map(|x| U192::from_dec_str(x).expect("error in divisor 128")).collect(),
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

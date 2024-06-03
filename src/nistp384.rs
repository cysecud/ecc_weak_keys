
pub mod scalar384;
pub mod field_p384;
pub mod implicit_group;

use self::scalar384::{U384,MathResult};
use self::field_p384::FieldP384;
use self::implicit_group::ImplicitP384;
use std::ops::Neg;

use crate::affine::AffinePoint;
use crate::FieldElement;
use crate::ellinit::{CurveParms, KeyPairs};
use crate::ellinit::Curve;
use crate::projective::ProjectivePoint;


#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
/*Structur of elliptic curve p192 */  
pub struct P384{
    pub q:U384,
    pub a:FieldP384,
    pub b:FieldP384,
}
impl <'a>CurveParms<'a> for P384 {
    const A:&'a str = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112316";
    const B:&'a str = "27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575";
    const XG:&'a str= "26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087";
    const YG:&'a str= "8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871";
}
impl Curve<FieldP384,U384,MathResult> for P384 {
    fn initialize()->Self {
        let q=U384::from_dec_str(FieldP384::PRIME).expect("error in Prime FieldP384 inizializing U384 curve");
        let a=FieldP384::new(U384::from_dec_str(Self::A).expect("Error in coefficient A"));
        let b=FieldP384::new(U384::from_dec_str(Self::B).expect("Error in coefficient B"));
            
        Self{a:a,b:b,q:q}
    }
    fn identity(&self)->AffinePoint<FieldP384>{
        AffinePoint{ 
            x: FieldP384::zero(), 
            y: FieldP384::one(), 
            infinity: 1 }
    }
    fn generator(&self)->ProjectivePoint<FieldP384> 
        {ProjectivePoint {
            x: FieldP384::new(U384::from_dec_str(Self::XG).expect("Error in x generator")),
            y: FieldP384::new(U384::from_dec_str(Self::YG).expect("Error in y generator")),
            z:FieldP384::one(),
            infinity: 0,}
    }
    fn ellnegate(&self,p:AffinePoint<FieldP384>)->AffinePoint<FieldP384> {
        AffinePoint { x: p.x, y: p.y.neg(), infinity: p.infinity }
    }
    fn ellnegate_proj(&self,p:ProjectivePoint<FieldP384>)->ProjectivePoint<FieldP384> {
        ProjectivePoint { x: p.x, y: p.y.neg(), z:p.z,infinity: p.infinity }
    }
     fn zn_prim_root_gen_order(&self)-> U384{
            U384::from_dec_str(ImplicitP384::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
     fn gen_order(&self)-> U384{
        U384::from_dec_str(ImplicitP384::PRIME).expect("error in generator's order P!")
    }

    fn ellisoncurve(&self,mut p:AffinePoint<FieldP384>)->bool {
        let a=FieldP384::new(U384::from_dec_str(Self::A).expect("error in A"));
        let b=FieldP384::new(U384::from_dec_str(Self::B).expect("error in B"));
        
        let x3=p.x.power(U384::from(3));
        let ax= a*p.x;
        let y2=x3+ax+b;
        let check:Option<U384>=U384::check_sqrt_mod_prime(y2.num, FieldP384::prime());
        if check.is_some() {
            if p.y.num==check.unwrap() || p.y.num==(FieldP384::prime()-check.unwrap()) {true} else{false}
        }else{false}
    }

    fn random (&self)->AffinePoint<FieldP384> {
        let mut x:FieldP384;
        let y:FieldP384;
        let a:FieldP384=FieldP384::new(U384::from_dec_str(Self::A).expect("error in A"));
        let b:FieldP384=FieldP384::new(U384::from_dec_str(Self::B).expect("error in B"));
        x=loop{
            x=FieldP384::rand_mod();
            let x3=FieldP384::power(&mut x,U384::from(3));
            let ax=a*x;
            let y2=x3+ax+b;
            let r:Option<U384>=U384::check_sqrt_mod_prime(y2.num,FieldP384::prime());
            if r.is_some() 
                {y=FieldP384::new(r.unwrap());
                    break x;}
        };
        AffinePoint::new(x,y)    }

    fn private_key(&self)->U384 {
        FieldP384::rand_mod().num    }

    fn key_pairs(&self)->KeyPairs<U384,AffinePoint<FieldP384>> {
            let sk=self.private_key();
            let pk = self.ellmul(&mut self.to_affine(self.generator()), sk);
            KeyPairs{sk,pk}
    }

    fn to_affine(&self,mut p:ProjectivePoint<FieldP384>)->AffinePoint<FieldP384> {
        let x_aff=p.x*FieldP384::inverse( &mut p.z);
        let y_aff=p.y*FieldP384::inverse(&mut p.z);
        AffinePoint::new(x_aff,y_aff)
    }

    fn jc_to_affine(&self,mut p:ProjectivePoint<FieldP384>)->AffinePoint<FieldP384> {
        let mut inv_z=p.z.inverse();
        let z=inv_z.square();
        let t=z*inv_z;
        let aff_x=p.x*z;
        let aff_y=p.y*t;
    AffinePoint::new(aff_x,aff_y)
    }

    fn proj_identity(&self)->ProjectivePoint<FieldP384> {
        ProjectivePoint{ 
            x: FieldP384::zero(), 
            y: FieldP384::one(),
            z: FieldP384::zero(), 
            infinity: 1 }
    }
    fn proj_elladd_jc(&self,mut p:ProjectivePoint<FieldP384>,mut q:ProjectivePoint<FieldP384>)->ProjectivePoint<FieldP384> {
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
    fn elladd(&self,aff_p:AffinePoint<FieldP384>,aff_q: AffinePoint<FieldP384>)->AffinePoint<FieldP384> {
            //let a4=FieldP384::new(U384::from_dec_str(Self::A).expect("error in A elladd"),aff_p.x.prime);
            //let q=aff_p.x.prime;
            if aff_p==Self::identity(&self) {return aff_q;}
            else if aff_q==Self::identity(self) {return aff_p;}
            else if aff_p==self.ellnegate(aff_q) {return  self.identity();}
    /*         else if q-U384::from(3)==a4.num {self.to_affine(self.proj_elladd_3(proj_p, proj_q))}
     */        else{
            let proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP384::one());
            let proj_q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldP384::one());
            self.jc_to_affine(self.proj_elladd_jc(proj_p, proj_q))}
    
}
fn proj_elldouble_jc(&self, proj_p:&mut ProjectivePoint<FieldP384>)->ProjectivePoint<FieldP384> {
    
    let a4:FieldP384=FieldP384::new(U384::from_dec_str(Self::A).expect("error in A proj_elldouble_jc"));

    let mut a=FieldP384::multiple(&(proj_p.x*(proj_p.y.square())),U384::from(4));
    let mut b:FieldP384;
    if FieldP384::prime()-U384::from(3)==a4.num {b=FieldP384::multiple(&((proj_p.x-proj_p.z.square())*(proj_p.x+proj_p.z.square())),U384::from(3));}
    else if a4.num==U384::zero(){b=FieldP384::multiple(&(proj_p.x.square()), U384::from(3))}
    else {b=FieldP384::multiple(&proj_p.x.square(),U384::from(3))+a4*proj_p.z.power(U384::from(4));}
    let x3=b.square()-a.double();
    let y3=-(proj_p.y.power(U384::from(4))).multiple(U384::from(8))+b*(a-x3);
    let z3=proj_p.y*proj_p.z.double();
    ProjectivePoint ::new(x3,y3, z3)
}
fn elldouble(&self,aff_p:&mut AffinePoint<FieldP384>)->AffinePoint<FieldP384> {
    let mut proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP384::one());
    self.jc_to_affine(self.proj_elldouble_jc(&mut proj_p))
}
fn ellmul(&self,p:&mut AffinePoint<FieldP384>,k:U384)->AffinePoint<FieldP384> {
    {

        /**************Montgomery ladder******************
[https://www.matthieurivain.com/files/jcen11b.pdf] 
Algorithm 3. */

let mut p0=self.proj_identity();
let mut p1= ProjectivePoint::new(p.x,p.y,FieldP384::one());
let bin=U384::to_binary(k);
for i in bin.iter(){
    if *i==U384::zero(){
        p1=self.proj_elladd_jc(p0, p1);

        p0=self.proj_elldouble_jc( &mut p0);}
    else {
        p0=self.proj_elladd_jc(p0, p1);
        p1=Self::proj_elldouble_jc(&self,&mut p1);}
    }
    self.jc_to_affine(p0)
}   

}
fn bsgs(&self,q:&mut AffinePoint<FieldP384>, d:U384)->MathResult {
        //This is the order of the generator point gen, it is calculated in pari-gp
	let ord:U384 = self.gen_order();
	//let mut z = FieldP384::znprimroot(&ord); /*generatore di Fp* dove p=ord(E,P)*/
	//This is a primitive root mod p, it is calculated in pari-gp
	let mut z= ImplicitP384::new(U384::from(self.zn_prim_root_gen_order()));
	let mut zd = ImplicitP384::power(&mut z,(ord-U384::one())/d);
	//calculate the inverse of the d-order generator
    let mut inv_zd=zd.inverse();
	let mut res:MathResult=None;
    let mut bs: Vec<(u128,AffinePoint<FieldP384>)>=Vec::new();
    let mut gs: Vec<(u128,AffinePoint<FieldP384>)>=Vec::new();
	let  m =U384::integer_sqrt(&d)+U384::one();//ceil(d)
	'outer: for i in 0..m.as_u128() {
        bs.push((i,self.ellmul(&mut self.jc_to_affine(self.generator()), ImplicitP384::power(&mut zd,U384::from(i)).num)));
        bs.sort_by_key(|key: &(u128, AffinePoint<FieldP384>)| key.1 );
		gs.push((i,self.ellmul(q,ImplicitP384::power(&mut inv_zd,m*U384::from(i)).num)));
        for j in 0..bs.len(){
			let search=bs.binary_search_by_key(&gs[j].1,|&(_a,b)|b);
            if search.is_ok()
			{
            let i =gs[j].0;
            let j=bs[search.unwrap()].0;
            println!("Baby Step Giant Step found a match:i={},j={}",i,j);
			let alpha=ImplicitP384::power(&mut zd,m*U384::from(i)+U384::from(j)%d);
			println!("The key is weak. Your secret key is {}",alpha.num);
            res=Some(alpha.num);
			break 'outer;
        }}}
			if res.is_some() {return res} else {return None;}	
}          

}


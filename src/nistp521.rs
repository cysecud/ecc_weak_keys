
pub mod scalar521;
pub mod field_p521;
pub mod implicit_group;

use self::scalar521::{U521,MathResult};
use self::field_p521::FieldP521;
use self::implicit_group::ImplicitP521;
use std::ops::Neg;

use crate::affine::AffinePoint;
use crate::FieldElement;
use crate::ellinit::CurveParms;
use crate::ellinit::Curve;
use crate::projective::ProjectivePoint;


#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
/*Structur of elliptic curve p192 */  
pub struct  P521{
    pub q:U521,
    pub a:FieldP521,
    pub b:FieldP521,
}

impl <'a>CurveParms<'a> for P521 {
    const A:&'a str = "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057148";
    const B:&'a str = "1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984";
    const XG:&'a str= "2661740802050217063228768716723360960729859168756973147706671368418802944996427808491545080627771902352094241225065558662157113545570916814161637315895999846";
    const YG:&'a str= "3757180025770020463545507224491183603594455134769762486694567779615544477440556316691234405012945539562144444537289428522585666729196580810124344277578376784";
    const P:&'a str = "6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449";
    const PRIME_ROOT:&'a str = "3";
}
impl Curve<FieldP521,U521,MathResult> for P521 {
    fn initialize()->Self {
        let q=U521::from_dec_str(FieldP521::PRIME).expect("error in Prime FieldP521 inizializing U521 curve");
        let a=FieldP521::new(U521::from_dec_str(Self::A).expect("Error in coefficient A"));
        let b=FieldP521::new(U521::from_dec_str(Self::B).expect("Error in coefficient B"));
            
        Self{a:a,b:b,q:q}
    }
    fn identity(&self)->AffinePoint<FieldP521>{
        AffinePoint{ 
            x: FieldP521::zero(), 
            y: FieldP521::one(), 
            infinity: 1 }
    }
    fn generator(&self)->ProjectivePoint<FieldP521> 
        {ProjectivePoint {
            x: FieldP521::new(U521::from_dec_str(Self::XG).expect("Error in x generator")),
            y: FieldP521::new(U521::from_dec_str(Self::YG).expect("Error in y generator")),
            z:FieldP521::one(),
            infinity: 0,}
    }
    fn ellnegate(&self,p:AffinePoint<FieldP521>)->AffinePoint<FieldP521> {
        AffinePoint { x: p.x, y: p.y.neg(), infinity: p.infinity }
    }
    fn ellnegate_proj(&self,p:ProjectivePoint<FieldP521>)->ProjectivePoint<FieldP521> {
        ProjectivePoint { x: p.x, y: p.y.neg(), z:p.z,infinity: p.infinity }
    }
     fn zn_prim_root_gen_order(&self)-> U521{
            U521::from_dec_str(Self::PRIME_ROOT).expect("error in primitive root of generator's order!")
    }
     fn gen_order(&self)-> U521{
        U521::from_dec_str(Self::P).expect("error in generator's order P!")
    }

    fn ellisoncurve(&self,mut p:AffinePoint<FieldP521>)->bool {
/*         let q:U521=U521::from_dec_str(Self::Q).expect("error in Q random");
 */     let a=FieldP521::new(U521::from_dec_str(Self::A).expect("error in A"));
        let b=FieldP521::new(U521::from_dec_str(Self::B).expect("error in B"));
        
        let x3=p.x.power(U521::from(3));
        let ax= a*p.x;
        let y2=x3+ax+b;
        let check:Option<U521>=U521::check_sqrt_mod_prime(y2.num, FieldP521::prime());
        if check.is_some() {
            if p.y.num==check.unwrap() || p.y.num==(FieldP521::prime()-check.unwrap()) {true} else{false}
        }else{false}
    }

    fn random (&self)->AffinePoint<FieldP521> {
        let mut x:FieldP521;
        let y:FieldP521;
        let a:FieldP521=FieldP521::new(U521::from_dec_str(Self::A).expect("error in A"));
        let b:FieldP521=FieldP521::new(U521::from_dec_str(Self::B).expect("error in B"));
        x=loop{
            x=FieldP521::rand_mod();
            let x3=FieldP521::power(&mut x,U521::from(3));
            let ax=a*x;
            let y2=x3+ax+b;
            let r:Option<U521>=U521::check_sqrt_mod_prime(y2.num,FieldP521::prime());
            if r.is_some() 
                {y=FieldP521::new(r.unwrap());
                    break x;}
        };
        AffinePoint::new(x,y)    }

    fn private_key(&self)->U521 {
        FieldP521::rand_mod().num    }

    fn public_key(&self,private_key:U521)->AffinePoint<FieldP521> {
        todo!()
    }

    fn to_affine(&self,mut p:ProjectivePoint<FieldP521>)->AffinePoint<FieldP521> {
        let x_aff=p.x*FieldP521::inverse( &mut p.z);
        let y_aff=p.y*FieldP521::inverse(&mut p.z);
        AffinePoint::new(x_aff,y_aff)
    }

    fn jc_to_affine(&self,mut p:ProjectivePoint<FieldP521>)->AffinePoint<FieldP521> {
        let mut inv_z=p.z.inverse();
        let z=inv_z.square();
        let t=z*inv_z;
        let aff_x=p.x*z;
        let aff_y=p.y*t;
    AffinePoint::new(aff_x,aff_y)
    }

    fn proj_identity(&self)->ProjectivePoint<FieldP521> {
        ProjectivePoint{ 
            x: FieldP521::zero(), 
            y: FieldP521::one(),
            z: FieldP521::zero(), 
            infinity: 1 }
    }
    fn proj_elladd_jc(&self,mut p:ProjectivePoint<FieldP521>,mut q:ProjectivePoint<FieldP521>)->ProjectivePoint<FieldP521> {
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
    fn elladd(&self,aff_p:AffinePoint<FieldP521>,aff_q: AffinePoint<FieldP521>)->AffinePoint<FieldP521> {
            //let a4=FieldP521::new(U521::from_dec_str(Self::A).expect("error in A elladd"),aff_p.x.prime);
            //let q=aff_p.x.prime;
            if aff_p==Self::identity(&self) {return aff_q;}
            else if aff_q==Self::identity(self) {return aff_p;}
            else if aff_p==self.ellnegate(aff_q) {return  self.identity();}
    /*         else if q-U521::from(3)==a4.num {self.to_affine(self.proj_elladd_3(proj_p, proj_q))}
     */        else{
            let proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP521::one());
            let proj_q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldP521::one());
            self.jc_to_affine(self.proj_elladd_jc(proj_p, proj_q))}
    
}
fn proj_elldouble_jc(&self, proj_p:&mut ProjectivePoint<FieldP521>)->ProjectivePoint<FieldP521> {
    
    let a4:FieldP521=FieldP521::new(U521::from_dec_str(Self::A).expect("error in A proj_elldouble_jc"));

    let mut a=FieldP521::multiple(&(proj_p.x*(proj_p.y.square())),U521::from(4));
    let mut b:FieldP521;
    if FieldP521::prime()-U521::from(3)==a4.num {b=FieldP521::multiple(&((proj_p.x-proj_p.z.square())*(proj_p.x+proj_p.z.square())),U521::from(3));}
    else if a4.num==U521::zero(){b=FieldP521::multiple(&(proj_p.x.square()), U521::from(3))}
    else {b=FieldP521::multiple(&proj_p.x.square(),U521::from(3))+a4*proj_p.z.power(U521::from(4));}
    let x3=b.square()-a.double();
    let y3=-(proj_p.y.power(U521::from(4))).multiple(U521::from(8))+b*(a-x3);
    let z3=proj_p.y*proj_p.z.double();
    ProjectivePoint ::new(x3,y3, z3)
}
fn elldouble(&self,aff_p:&mut AffinePoint<FieldP521>)->AffinePoint<FieldP521> {
    let mut proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldP521::one());
    self.jc_to_affine(self.proj_elldouble_jc(&mut proj_p))
}
fn ellmul(&self,p:&mut AffinePoint<FieldP521>,k:U521)->AffinePoint<FieldP521> {
    {

        /**************Montgomery ladder******************
[https://www.matthieurivain.com/files/jcen11b.pdf] 
Algorithm 3. */

let mut p0=self.proj_identity();
let mut p1= ProjectivePoint::new(p.x,p.y,FieldP521::one());
let bin=U521::to_binary(k);
for i in bin.iter(){
    if *i==U521::zero(){
        p1=self.proj_elladd_jc(p0, p1);

        p0=self.proj_elldouble_jc( &mut p0);}
    else {
        p0=self.proj_elladd_jc(p0, p1);
        p1=Self::proj_elldouble_jc(&self,&mut p1);}
    }
    self.jc_to_affine(p0)
}   

}
fn bsgs(&self, q:&mut AffinePoint<FieldP521>, d:U521)->MathResult {
        //This is the order of the generator point gen, it is calculated in pari-gp
	let ord:U521 = self.gen_order();
	//let mut z = FieldP521::znprimroot(&ord); /*generatore di Fp* dove p=ord(E,P)*/
	//This is a primitive root mod p, it is calculated in pari-gp
	let mut z= ImplicitP521::new(U521::from(self.zn_prim_root_gen_order()));
	let mut zd = ImplicitP521::power(&mut z,(ord-U521::one())/d);
	//calculate the inverse of the d-order generator
    let mut inv_zd=zd.inverse();
	let mut res:MathResult=None;
    let mut bs: Vec<(u128,AffinePoint<FieldP521>)>=Vec::new();
    let mut gs: Vec<(u128,AffinePoint<FieldP521>)>=Vec::new();
	let  m =U521::integer_sqrt(&d)+U521::one();//ceil(d)
	'outer: for i in 0..m.as_u128() {
        bs.push((i,self.ellmul(&mut self.jc_to_affine(self.generator()), ImplicitP521::power(&mut zd,U521::from(i)).num)));
        bs.sort_by_key(|key: &(u128, AffinePoint<FieldP521>)| key.1 );
		gs.push((i,self.ellmul(q,ImplicitP521::power(&mut inv_zd,m*U521::from(i)).num)));
        for j in 0..bs.len(){
			let search=bs.binary_search_by_key(&gs[j].1,|&(_a,b)|b);
            if search.is_ok()
			{
               println!("bs is {:?}",bs[search.unwrap()]);
                println!("gs is {:?}",gs[j]);
            let i =gs[j].0;
            let j=bs[search.unwrap()].0;
            println!("Baby Step Giant Step found a match:i={},j={}",i,j);
			let alpha=ImplicitP521::power(&mut zd,m*U521::from(i)+U521::from(j)%d);
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

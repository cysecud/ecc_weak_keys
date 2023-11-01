use crate::projective::ProjectivePoint;
use crate::scalar::{U256,MathResult};
use crate::field::FieldPoint;
use crate::affine::AffinePoint;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
pub struct Ellinit{
    pub a:FieldPoint,
    pub b:FieldPoint,
    pub q:U256
}
pub type Ell=Ellinit;
pub trait CurveParms<'a>{    
    //required//
    
    const Q:&'a str;
    const A:&'a str;
    const B:&'a str;
    fn generator(&self)->ProjectivePoint;
    fn identity(&self)->AffinePoint;
    fn zn_prim_root_gen_order(&self)->U256;
    fn gen_order(&self)->U256;

    fn ellisoncurve(&self,mut p:AffinePoint)->bool {
        let q:U256=U256::from_dec_str(Self::Q).expect("error in Q random");
        let a:FieldPoint=FieldPoint::new(U256::from_dec_str(Self::A).expect("error in Q random"),q);
        let b:FieldPoint=FieldPoint::new(U256::from_dec_str(Self::B).expect("error in Q random"),q);
        
        let x3=FieldPoint::power( &mut p.x, U256::from(3));
        let ax= a*p.x;
        let y2=x3+ax+b;
        let check:Option<U256>=U256::check_sqrt_mod_prime(y2.num, q);
        if check.is_some() {
            if p.y.num==check.unwrap() || p.y.num==(p.y.prime-check.unwrap()) {true} else{false}
        }else{false}
    }

    fn random (&self)->AffinePoint{
        let mut x:FieldPoint;
        let y:FieldPoint;
        let q:U256=U256::from_dec_str(Self::Q).expect("error in Q random");
        let a:FieldPoint=FieldPoint::new(U256::from_dec_str(Self::A).expect("error in Q random"),q);
        let b:FieldPoint=FieldPoint::new(U256::from_dec_str(Self::B).expect("error in Q random"),q);

        x=loop{
            x=FieldPoint::rand_mod(q);
            let x3=FieldPoint::power(&mut x,U256::from(3));
            let ax=a*x;
            let y2=x3+ax+b;
            let r:Option<U256>=U256::check_sqrt_mod_prime(y2.num,x.prime);
            if r.is_some() 
                {y=FieldPoint::new(r.unwrap(), x.prime);
                    break x;}
        };
        AffinePoint::new(x,y)
    }
    fn private_key(&self)->U256{
        let q:U256=U256::from_dec_str(Self::Q).expect("error in Q random");
        FieldPoint::rand_mod(q).num
    }
    fn public_key(&self,private_key:U256)->AffinePoint{
        self.ellmul(&mut self.to_affine(self.generator()), private_key)
    }
    fn to_affine(&self,mut p:ProjectivePoint)->AffinePoint{
        let x_aff=p.x*FieldPoint::inverse( &mut p.z);
        let y_aff=p.y*FieldPoint::inverse(&mut p.z);
        AffinePoint::new(x_aff,y_aff)
    }
    fn jc_to_affine(&self,mut p:ProjectivePoint)->AffinePoint{
        let mut inv_z=p.z.inverse();
        let z=inv_z.square();
        let t=z*inv_z;
        let aff_x=p.x*z;
        let aff_y=p.y*t;
    AffinePoint::new(aff_x,aff_y)
}
    fn proj_identity(&self)->ProjectivePoint {
        ProjectivePoint{ 
            x: FieldPoint{num:U256::zero(),prime:self.identity().x.prime}, 
            y: FieldPoint{num:U256::zero(),prime:self.identity().x.prime},
            z: FieldPoint{num:U256::zero(),prime:self.identity().x.prime}, 
            infinity: 1 }
    }
    fn elladd(&self,aff_p:AffinePoint,aff_q: AffinePoint)->AffinePoint{
        
        //let a4=FieldPoint::new(U256::from_dec_str(Self::A).expect("error in A elladd"),aff_p.x.prime);
        //let q=aff_p.x.prime;

        if aff_p==Self::identity(&self) {return aff_q;}
        else if aff_q==Self::identity(self) {return aff_p;}
        else if aff_p==aff_q.negate() {return  self.identity();}
/*         else if q-U256::from(3)==a4.num {self.to_affine(self.proj_elladd_3(proj_p, proj_q))}
 */        else{
        let proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldPoint::new(U256::from(1),U256::from_dec_str(&Self::Q).expect("error in Q")));
        let proj_q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldPoint::new(U256::from(1),U256::from_dec_str(&Self::Q).expect("error in Q")));
        self.jc_to_affine(self.proj_elladd_jc(proj_p, proj_q))}
    
}

fn elldouble(&self,aff_p:&mut AffinePoint)->AffinePoint{
        let mut proj_p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldPoint::new(U256::from(1),U256::from_dec_str(&Self::Q).expect("error in Q")));
        self.jc_to_affine(self.proj_elldouble_jc(&mut proj_p))
}
fn proj_elldouble_jc(&self, proj_p:&mut ProjectivePoint)->ProjectivePoint{

        let a4:FieldPoint=FieldPoint::new(U256::from_dec_str(Self::A).expect("error in A proj_elldouble_jc"),proj_p.x.prime);
        let q=U256::from_dec_str(Self::Q).expect("Error in Q in projective proj_elldouble_jc");

        let a=FieldPoint::multiple(&(proj_p.x*(proj_p.y.square())),U256::from(4));
        let mut b:FieldPoint;
        if q-U256::from(3)==a4.num {b=FieldPoint::multiple(&((proj_p.x-proj_p.z.square())*(proj_p.x+proj_p.z.square())),U256::from(3));}
        else if a4.num==U256::zero(){b=FieldPoint::multiple(&(proj_p.x.square()), U256::from(3))}
        else {b=FieldPoint::multiple(&proj_p.x.square(),U256::from(3))+a4*proj_p.z.power(U256::from(4));}
        let x3=a.double().negate()+b.square();
        let y3=FieldPoint::multiple(&(proj_p.y.power(U256::from(4))), U256::from(8)).negate()+b*(a-x3);
        let z3=proj_p.y*proj_p.z.double();
        ProjectivePoint ::new(x3,y3, z3)

}
fn proj_elladd_jc(&self,mut p:ProjectivePoint,mut q:ProjectivePoint)->ProjectivePoint{
    /*Jacobian and Chudnovsky Jacobian coordinates
    With Jacobian coordinates the curve E is given by
                Y^2 = X^3 + a4XZ 4 + a6Z^6 .
    The point (X1 : Y1 : Z1 ) on E corresponds to the affine point (X1 /Z1^2 , Y1 /Z1^3 ) when Z1 is not 0 and to
    the point at infinity otherwise. The opposite of (X1 : Y1 : Z1 ) is (X1 : −Y1 : Z1 ). 
    The complexities are 12M + 4S for an addition and 4M + 6S for a doubling. If one of the points is
    given in the form (X1 : Y1 : 1) the costs for addition reduce to 8M + 3S.*/
    if p==self.proj_identity() {return q;} 
    else if q==self.proj_identity() {return p;}
    else if p==q.negate() {return self.proj_identity();}
    else{
        let a=p.x*q.z.square();
        let b=q.x*p.z.square();
        let c=p.y*q.z.square()*q.z;
        let d=q.y*p.z.square()*p.z;
        let mut e=b-a;
        let mut f=d-c;
        let x3=(e.square()*e).negate()+(a*e.square().double()).negate()+f.square();
        let y3=(c*e.square()*e).negate()+f*(a*e.square()-x3);
        let z3=p.z*q.z*e;
        ProjectivePoint::new(x3,y3,z3)
    }
}
fn ellmul(&self,p:&mut AffinePoint,k:U256)->AffinePoint{

            /**************Montgomery ladder******************
    [https://www.matthieurivain.com/files/jcen11b.pdf] 
    Algorithm 3. */

    let mut p0=ProjectivePoint{x:Self::identity(&self).x,y:Self::identity(&self).y,z:Self::identity(&self).x,infinity:1};
    let mut p1= ProjectivePoint::new(p.x,p.y,FieldPoint::new(U256::one(),p.x.prime));
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
        /*Implicit baby-step-giant-step algorithm */
fn bsgs(&self, gen:&mut AffinePoint,q:&mut AffinePoint, d:U256)->MathResult {
    //This is the order of the generator point gen, it is calculated in pari-gp
	let ord:U256 = self.gen_order();
	//let mut z = FieldPoint::znprimroot(&ord); /*generatore di Fp* dove p=ord(E,P)*/
	//This is a primitive root mod p, it is calculated in pari-gp
	let mut z= FieldPoint::new(U256::from(self.zn_prim_root_gen_order()),ord);
	let mut zd = FieldPoint::power(&mut z,(ord-U256::one())/d);
	//calculate the inverse of the d-order generator
    let mut inv_zd=zd.inverse();
	let mut res:MathResult=None;
    let mut bs: Vec<(u128,AffinePoint)>=Vec::new();
    let mut gs: Vec<(u128,AffinePoint)>=Vec::new();
	let  m =U256::integer_sqrt(&d)+U256::one();//ceil(d)
	'outer: for i in 0..m.as_u128() {
        bs.push((i,self.ellmul(gen, FieldPoint::power(&mut zd,U256::from(i)).num)));
        bs.sort_by_key(|key: &(u128, AffinePoint)| key.1 );
		gs.push((i,self.ellmul(q,FieldPoint::power(&mut inv_zd,m*U256::from(i)).num)));
        for j in 0..bs.len(){
			let search=bs.binary_search_by_key(&gs[j].1,|&(_a,b)|b);
            if search.is_ok()
			{
               println!("bs is {:?}",bs[search.unwrap()]);
                println!("gs is {:?}",gs[j]);
            let i =gs[j].0;
            let j=bs[search.unwrap()].0;
            println!("Baby Step Giant Step found a match:i={},j={}",i,j);
			let alpha=FieldPoint::power(&mut zd,m*U256::from(i)+U256::from(j)%d);
			println!("The key is weak. Your secret key is {}",alpha.num);
            res=Some(alpha.num);
			break 'outer;
        }}}
			if res.is_some() {return res} else {return None;}	
}          

}


impl Ellinit {
    pub fn new(a:FieldPoint,b:FieldPoint,q:U256)-> Self {
        Ellinit { a: a, b: b, q: q}
    }
    }




pub(crate) use core::ops::{Add,Neg};
use crate::ellinit::Ellinit;
use crate::scalar::U256;
use crate::field::FieldPoint;

/*this structur in in the lib.rs
It contains the parameters of the elliptic curve p224 */
/*----------------------Projective point structure------------*/
    
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
pub struct ProjectivePoint {
    pub x: FieldPoint,
    pub y: FieldPoint,
    pub z: FieldPoint,

    /// Is this point the point at infinity? 0 = no, 1 = yes
    ///
    pub infinity: u8,
}
impl ProjectivePoint {
    /// Create a new [`AffinePoint`] with the given coordinates.
    pub fn new(x: FieldPoint, y: FieldPoint,z:FieldPoint) -> Self {
        Self { x, y, z,infinity: 0 }
    }  
    pub fn negate(&self) -> Self {
        ProjectivePoint  {
            x: self.x,
            y: self.y.negate(),
            z:self.z,
            infinity: self.infinity,
        }
    }  
}
impl Neg for ProjectivePoint {
    type Output = Self;

    fn neg(self) -> Self::Output {
        ProjectivePoint {
            x: self.x,
            y: self.y.negate(),
            z: self.z,
            infinity: self.infinity,
        }
    }
}

/*----------------Arithmetic for projective points---------*/
/* impl Add for ProjectivePoint {
    type Output = Self;
    fn add( self, other: Self) -> Self {	
	if self==ProjectivePoint::identity() {return other;}
	else if other==ProjectivePoint::identity() {return self;}
    else if self == other.neg() {return ProjectivePoint::identity()}
        else { 
		
		let mut a= other.y*self.z+(self.y*other.z).negate(); 
		let mut b= other.x*self.z+(self.x*other.z).negate();
		let c= FieldPoint::power(&mut a,U256::from(2u128))*self.z*other.z+FieldPoint::power(&mut b,U256::from(3u128)).negate()+(FieldPoint::power(&mut b,U256::from(2u128))*self.x*other.z).negate()+(FieldPoint::power(&mut b,U256::from(2u128))*self.x*other.z).negate();
		let x3=b*c;
		let y3=a*(FieldPoint::power(&mut b,U256::from(2u128))*self.x*other.z+c.negate())+(FieldPoint::power(&mut b,U256::from(3u128))*self.y*other.z).negate();
		let z3= FieldPoint::power(&mut b,U256::from(3u128))*self.z*other.z;
		
         /*This is an implementation for addition that works if coefficient a of the elliptic curve is a=-3 */
  //       let p=ProjectivePoint::new(aff_p.x, aff_p.y, FieldPoint::new(U256::from(1),self.q));
//let q=ProjectivePoint::new(aff_q.x, aff_q.y, FieldPoint::new(U256::from(1),self.q));
/* let a6=FieldPoint::new(U256::from_dec_str(NISTp224::B).expect("error in B coefficient in projective addition formula"),U256::from_dec_str(NISTp224::Q).expect("Error in prime Q in projective addition formula"));

    let t0: FieldPoint= self.x*other.x;    //1
    let t1: FieldPoint =self.y*other.y;    //2
    let t2: FieldPoint =self.z*other.z;    //3
    let t3: FieldPoint=self.x+self.y;     //4
    let t4: FieldPoint=other.x+other.y;    //5
    let t3: FieldPoint=t3*t4;              //6
    let t4: FieldPoint= t0 + t1;           //7
    let t3: FieldPoint=t3-t4; //8
    let t4: FieldPoint= self.y+self.z; //9
    let x3: FieldPoint = other.y+other.z; //10
    let t4: FieldPoint = t4*x3;//11
    let x3: FieldPoint = t1 + t2; //12
    let t4: FieldPoint = t4-x3; //13
    let x3: FieldPoint = self.x + self.z; //14
    let y3: FieldPoint =other.x + other.z;//15
    let x3: FieldPoint = x3*y3; //16
    let y3: FieldPoint = t0 + t2; //17
    let y3: FieldPoint = x3 -y3;//18
    let z3: FieldPoint = a6*t2; //19
    let x3: FieldPoint = y3-z3; //20 
    let z3: FieldPoint = x3 + x3; //21
    let x3: FieldPoint=x3 + z3; //22
    let z3: FieldPoint=t1 - x3; //23
    let x3: FieldPoint = t1 + x3; //24
    let y3: FieldPoint= a6*y3;   //25               //25
    let t1: FieldPoint= t2 + t2;   //26
    let t2: FieldPoint = t1 + t2;  //27
    let y3: FieldPoint =y3-t2; //28
    let y3: FieldPoint =y3-t0; //29
    let t1: FieldPoint = y3 + y3; //30
    let y3: FieldPoint = t1+y3; //31
    let t1: FieldPoint = t0+t0; //32
    let t0: FieldPoint =t1 + t0; //33
    let t0: FieldPoint = t0-t2; //34
    let t1: FieldPoint = t4*y3;//35
    let t2: FieldPoint = t0*y3;//36
    let y3: FieldPoint =x3* z3;//37
    let y3: FieldPoint =y3+t2;//38
    let x3: FieldPoint = t3*x3;//39
    let x3: FieldPoint = x3-t1;//40
    let z3: FieldPoint = t4*z3;//41
    let t1: FieldPoint = t3*t0;//42
    let z3: FieldPoint = z3+t1;//43
   
 */   ProjectivePoint::new(x3,y3,z3)
        }    
    }
} */

impl ProjectivePoint {
    ///Double point
    pub fn double(&mut self,ell:Ellinit)->Self {
    //Here we need the coefficients of the curve y^2z=x^3+axz^2+bz^3.
    //This works...
    let a4 = ell.a;
	
    let mut a =a4*FieldPoint::power( &mut self.z,U256::from(2))+FieldPoint::multiple(&FieldPoint::power(&mut self.x,U256::from(2)),U256::from(3));
    let mut b = self.y*self.z;
    let mut c = b*self.x*self.y;
    let d = FieldPoint::power(&mut a,U256::from(2))+FieldPoint::multiple(&c,U256::from(8)).negate();
    	let x3= FieldPoint::multiple(&mut (b*d),U256::from(2));
    	let y3=	a*(FieldPoint::multiple(&mut c,U256::from(4))+d.negate())+FieldPoint::multiple(&FieldPoint::power(&mut (self.y*b),U256::from(2)),U256::from(8)).negate();
    	let z3=FieldPoint::multiple(&FieldPoint::power(&mut b,U256::from(3)),U256::from(8));
    	
        ProjectivePoint::new(x3,y3,z3)
    /*Implements point doubling for curves with `a = -3`
    Implements the exception-free point doubling formula from [Renes-Costello-Batina 2015]
    (Algorithm 6). The comments after each line indicate which algorithm
    steps are being performed. */
   /*  if self==ProjectivePoint::identity(){return  ProjectivePoint::identity();}
    else {
    let a6=FieldPoint::new(U256::from_dec_str(NISTp224::B).expect("error in B coefficient in projective addition formula"),U256::from_dec_str(NISTp224::Q).expect("Error in prime Q in projective addition formula"));
    let mut p=ProjectivePoint::new(self.x, self.y, FieldPoint::new(U256::from(1),U256::from_dec_str(NISTp224::Q).expect("Error in prime Q in projective addition formula")));

    let t0 = FieldPoint::square(&mut p.x); //1
    let t1 = FieldPoint::square(&mut p.y);//2 
    let t2=FieldPoint::square(&mut p.z); //3
    let t3 = p.x * p.y;//4
    let t3 = t3+t3; //5
    let z3 = p.x*p.z; //6
    let z3 = z3+z3;//7
    let y3 = a6*t2;//8
    let y3 = y3 - z3; //9
    let x3 = y3 + y3;
    let y3 = x3+ y3;
    let x3 = t1 - y3;
    let y3 = t1 + y3;
    let y3 = x3* y3;
    let x3 = x3*t3;
    let t3 = t2 + t2;
    let t2 = t2 + t3;
    let z3 = a6*z3;
    let z3 = z3-t2;
    let z3 = z3 -t0;
    let t3 = z3 + z3;
    let z3 = z3 + t3; 
    let t3 = t0 + t0;
    let t0 = t3 + t0;
    let t0 = t0 -t2;
    let t0 = t0*z3;
    let y3 = y3 + t0;
    let t0: FieldPoint = p.y*p.z;
    let t0 = t0 + t0;
    let z3 = t0 *z3;
    let x3=x3-z3;
    let z3 = t0*t1;
    let z3 = z3 + z3;
    let z3 = z3 + z3;
ProjectivePoint::new(x3,y3,z3)}
 */
   }
   /* pub fn multiple(ell:Ellinit,aff_p:&mut AffinePoint, k:U256) -> Self {        
        /**************Montgomery ladder******************
    [https://www.matthieurivain.com/files/jcen11b.pdf] 
    Algorithm 3. */
    let mut p0=ProjectivePoint::identity();
    let mut p1= ProjectivePoint::new(aff_p.x,aff_p.y,FieldPoint::new(U256::one(),aff_p.x.prime));
    let bin=U256::to_binary(k);
    for i in bin.iter(){
        if *i==U256::zero(){
            p1=p1+p0;

            p0=p0.double(ell);

}
        else {
            p0=p0+p1;
            p1=p1.double(ell);

}

        }
        p0
    }
 */
/* pub fn to_affine(self) ->AffinePoint{
    if self==ProjectivePoint::identity(){AffinePoint::identity()}
    else{
    let x_aff=self.x*FieldPoint::inverse( &self.z);
    let y_aff=self.y*FieldPoint::inverse(&self.z);
    AffinePoint::new(x_aff,y_aff)}
} */
}


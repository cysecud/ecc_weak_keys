use elliptic_curves::ellinit::{Curve,CurveParms};
use elliptic_curves::FieldElement;
use elliptic_curves::nistp192::{field_p192::FieldP192,scalar192::U192,P192};
use elliptic_curves::nistp224::{field_p224::FieldP224,scalar224::U224,P224};
use elliptic_curves::nistp384::{field_p384::FieldP384,scalar384::U384,P384};
use elliptic_curves::nistp521::{field_p521::FieldP521,scalar521::U521,P521};

fn main () {   
        /*Test Field192 */ 
let mut a=FieldP192::rand_mod();
let  c=FieldP192::rand_mod();
let sum=a*c;
let inv =a.inverse();
println!("inverse of {} is {}",a.num,inv.num);
println!("Field192 prime base is {}",FieldP192::prime());
println!("{}*{}={}",a.num,c.num,sum.num);
let p192=P192::initialize();
println!("gen p192 {:?}",p192.generator());
println!("identity p192 {:?}",p192.identity());
let mut x=p192.random();
let mut ison:bool=p192.ellisoncurve(x);
    if ison {println!("Point is on curve")} else {println!("Point NOT on the curve")};
let y=p192.random();
let sum = p192.elladd(x, y);
println!("{:?}+{:?}={:?}",x,y,sum);
let doub=p192.elldouble(&mut x);
println!("double of {:?}={:?}",x,doub);
let k=U192::from_dec_str("35692188470218831599429521488958156389612077476752115046").expect("error");
let mut mult=p192.ellmul(&mut p192.to_affine(p192.generator()), k);
println!("Point to be tested is {:?}",mult);
let d=U192::from(16);
let alpha=p192.bsgs(&mut p192.to_affine(p192.generator()), &mut mult, d);

    /*Test Field224 */ 
let mut a=FieldP224::rand_mod();
let  c=FieldP224::rand_mod();
let sum=a*c;
let inv =a.inverse();
println!("inverse of {} is {}",a.num,inv.num);
println!("Field192 prime base is {}",FieldP224::prime());
println!("{}*{}={}",a.num,c.num,sum.num);
let p224=P224::initialize();
println!("gen p192 {:?}",p224.generator());
println!("identity p192 {:?}",p224.identity());
let mut x=p224.random();
let mut ison:bool=p224.ellisoncurve(x);
    if ison {println!("Point is on curve")} else {println!("Point NOT on the curve")};
let y=p224.random();
let sum = p224.elladd(x, y);
println!("{:?}+{:?}={:?}",x,y,sum);
let doub=p224.elldouble(&mut x);
println!("double of {:?}={:?}",x,doub);
let k=U224::from_dec_str("23430269157991683102913750948896682204368511013412913881703591315224").expect("error");
let mut mult=p224.ellmul(&mut p224.to_affine(p224.generator()), k);
println!("Point to be tested is {:?}",mult);
let d=U224::from(12);
let alpha=p224.bsgs(&mut p224.to_affine(p224.generator()), &mut mult, d);
/*          
            /*Test Field256k1 */ 
let mut a=FieldP256k1::rand_mod();
let c=FieldP256k1::rand_mod();
let sum=a*c;
let inv =a.inverse();
println!("inverse of {} is {}",a.num,inv.num);
println!("Field192 prime base is {}",FieldP256k1::prime());
println!("{}*{}={}",a.num,c.num,sum.num);
    /*Test Field384 */ 
let mut a=FieldP384::rand_mod();
let c=FieldP384::rand_mod();
let sum=a*c;
let inv =a.inverse();
println!("inverse of {} is {}",a.num,inv.num);
println!("Field192 prime base is {}",FieldP384::prime());
println!("{}*{}={}",a.num,c.num,sum.num);
    /*Test Field521  */ 
let mut a=FieldP521::rand_mod();
let c=FieldP521::rand_mod();
let sum=a*c;
let inv =a.inverse();
println!("inverse of {} is {}",a.num,inv.num);
println!("Field192 prime base is {}",FieldP521::prime());
println!("{}*{}={}",a.num,c.num,sum.num);
    
 */
/* 
println!("{} mod {} is {}",rand,prime,rand%prime);
println!("{} mod {} is {}",rand,prime,rand.mod_fast(prime));
let mut x =FieldPoint::rand_mod();
let mut y=FieldPoint::rand_mod();
    println!("field point x = {}",x.num);
    println!("field point x = {}",y.num);
    let sum =x+y;
    let mul=x*y;
    let pow=x.power(y.num);
    println!("sum is {:?}",sum);
    println!("mul is {:?}",mul);
    println!("pow is {:?}",pow);

 */
    
 /*    
println!("quadratic residue mod 101");
    let p:U256= U256::from(101);
for i in 1..(p-1).as_usize() 
{if U256::quad_res(U256::from(i as u8), p){println!("{}",i);} }
println!("legendre mod 101");
for i in 1..(p-1).as_usize() {println!("legendre of {} with respect to {} is {}", i,p,U256::legendre(i.into(),p))}
 */

}

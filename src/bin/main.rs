use elliptic_curves::affine::AffinePoint;
use elliptic_curves::ellinit::CurveParms;
use elliptic_curves::nistp192::P192;
use elliptic_curves::nistp224::P224;
use elliptic_curves::nistp256k1::P256k1;
use elliptic_curves::nistp256r1::P256r1;
use elliptic_curves::nistsm2::Sm2;
use elliptic_curves::scalar::U256;
use elliptic_curves::field::FieldPoint;
fn main (){
    
    /*INITIALISATION OF THE CURVE P256K1 */
    let p_256k1= P256k1::initialize();
    println!("Curve P256K1 has been initialized!");
    let p_224= P224::initialize();
    println!("Curve p224 has been initialized!");
    let p_192= P192::initialize();
    println!("Curve p192 has been initialized!");
    let p_256r1= P256r1::initialize();
    println!("Curve p256r1 has been initialized!");
    let sm2= Sm2::initialize();
    println!("Curve Sm2 has been initialized!");

    let  point1=p_256k1.random();
    println!("point1 is {:?}",point1);
    let mut ison:bool=p_256k1.ellisoncurve(point1);
    if ison {println!("Point1 is on curve")} else {println!("Point1 NOT on the curve")};

    /* let  point2=p_224.random();
    println!("point2 is {:?}",point2);
    ison=p_224.ellisoncurve(point2);
    if ison {println!("Point2 is on curve")} else {println!("Point2 NOT on the curve")};

    let  point3=p_192.random();
    println!("point3 is {:?}",point3);
    ison=p_192.ellisoncurve(point3);
    if ison {println!("Point3 is on curve")} else {println!("Point3 NOT on the curve")};
    
    let  point4=p_256r1.random();
    println!("point4 is {:?}",point4);
    ison=p_256r1.ellisoncurve(point4);
    if ison {println!("Point4 is on curve")} else {println!("Point4 NOT on the curve")};

    let  point5=sm2.random();
    println!("point4 is {:?}",point5);
    ison=sm2.ellisoncurve(point5);
    if ison {println!("Point5 is on curve")} else {println!("Point5 NOT on the curve")};
 */
    let k=U256::from_dec_str("39283872389482783983661720076578799332063276928947188074161284086913226340463").expect("error!!!");
    let mut mul1=p_256k1.ellmul(&mut p_256k1.to_affine(p_256k1.generator()),k);
/*     let mut mul2=p_256r1.ellmul(&mut p_256r1.to_affine(p_256r1.generator()),k);
    let mut mul3=p_224.ellmul(&mut p_224.to_affine(p_224.generator()),k);
    let mut mul3=p_192.ellmul(&mut p_192.to_affine(p_192.generator()),k);
    let mut mul4=sm2.ellmul(&mut sm2.to_affine(sm2.generator()),k);

    println!("multiple of point 1 is {:?}",mul1);
    println!("multiple of point 2 is {:?}",mul2);
    println!("multiple of point 3 is {:?}",mul3);
    println!("multiple of point 4 is {:?}",mul4);
 */
    //let pk=p_256k1.private_key();
   // println!("A random generated private key in P256K1 is {}", pk);

    let mut gen:AffinePoint=p_256k1.to_affine(p_256k1.generator());
    println!("The generator of the curve P256K1 is {:?}",gen);
    //TEST OF WEAKNESS ON THE CURVE P256K1//
    /* let xq=FieldPoint::new(U256::from_dec_str("100760202697161893004335214126591116800117319792545458764085267675326325395621").expect("error"),p_256k1.q);
    let yq=FieldPoint::new(U256::from_dec_str("75193444318165031146359304621062797862272142296678797285916994295833810377664").expect("error"),p_256k1.q);
     */

    /* let xq=FieldPoint::new(U256::from_dec_str("88294711626906733189146820098032765429493710408815299173986011815090993120712").expect("error"),p_256k1.q);
    let yq=FieldPoint::new(U256::from_dec_str("43570037339456764009332097386535472658659524190102285899685553776705216363095").expect("error"),p_256k1.q);
    
    let mut public_key:AffinePoint=AffinePoint { x: xq, y: yq, infinity: 0 };
    println!("public key to be tested {:?}",public_key);
    let ison:bool=p_256k1.ellisoncurve(public_key);
    if ison {println!("Point is on curve")} else {println!("Point NOScalar on the curve")};
    */
    let div:U256 = U256::from(15144u128);
    println!("The divisor used for the test is {:?}",div);
    p_256k1.bsgs(&mut gen, &mut mul1, div);
    /*INIZIALIE THE CURVE P192 */  
    }
    
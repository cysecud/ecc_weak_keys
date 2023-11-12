use crate::Scalar;
use rand::Rng;
use uint::construct_uint;
construct_uint!{
    pub struct U384(12);}
pub type MathResult = Option<U384>;
type RootsResult = [Option<U384>;2];
impl <'a>Scalar<'a>for U384{
    const ORDER:&'a str = "384";
    const SIZE:usize = 12;
}

impl U384 {
    /*Algorithm 1 Crandall [https://eprint.iacr.org/2022/411.pdf] 
    Pseudo Mersenne Numbers are integers of the form q = 2^l − c,
    where c is ”small”. An algorithm to computhe their modular reduction is introduced by Crandall*/

    pub fn mod_fast(mut self, rhs: Self) -> Self{
        let l=Self::from_dec_str(Self::ORDER).expect("error fast modular");
        let p=Self::from(2).pow(l);
        let c=p-rhs;
        while self>=Self::from(2)*rhs {
            let a0=self&(p-1);
            let a1=self>>l;
            self=c*a1+a0;}
            if self<rhs{return self} else {return self-rhs}
 }
    pub fn random()->Self {
        let mut rng = rand::thread_rng();
        let mut v_rand = [0u64; Self::SIZE];
        rng.fill(&mut v_rand);
        Self(v_rand)
    }
    pub fn to_binary(mut self)->Vec<Self>{
        let mut bin: Vec<Self>=Vec::new();
        while self!=Self::zero() {
            bin.push(self%Self::from(2));
            self=self/Self::from(2);};
        bin= bin.into_iter().rev().collect();
       bin}
    pub fn bin_to_dec(bin:Vec<Self>)->Self{
        let mut dec:Self=Self::zero();
        for i in 0..bin.len(){
            dec=dec+bin[i]*Self::pow(Self::from(2), (bin.len()-1-i).into());
           
        }
        dec
    }
    pub const TAB2:[i8;8]=[0,1,0,-1,0,-1,0,1];
    pub fn kronecker (mut a:Self,mut b:Self)->i8 {
        /*inizialize k in this scope */
        let mut k:i8;

        /*Test b = 0 */
    
if b==Self::zero() {if a!=Self::one() {return 0i8;} else if  a==Self::one() {return 1i8;} }

/*Remove 2 from b: v will be the even part of b */
if a%2==Self::zero() && b%2==Self::zero() {return 0i8;} 
        else {
            /*initialize v = 0 */
            let mut v:Self=Self::zero();
            /*calulate the even part of b */
            while b%2==Self::zero() {
                v+=Self::one();
                b=b/2;                
            }
            if v%2==Self::zero(){k=1;} else {k=Self::TAB2[(a&Self::from(7)).as_usize()];}
        
        }
        /*now b must me odd */
/*         println!("b is odd? {}",b); 
 */    /*Recuce size once */
        a=a%(b);
    /*Finished? Loop.. */
    loop{
        if a==Self::zero() {if b>Self::one() {break 0i8;} else if b==Self::one() {break k;}}
    /*Remove power of 2 from a */
    /*initialize v = 0 */
            let mut v=Self::zero();
            /*calulate the even part of b */
            while a%2==Self::zero() {
                v+=Self::one();
                a=a/Self::from(2);                
            }
        if v%2==Self::one() {k=Self::TAB2[(b&Self::from(7)).as_usize()]*k;}
        /*subtract and apply reciprocity */
        /* a and b are both odd now */
        if b>a {let r:Self=b-a; 
                if (a&b&Self::from(2))!=Self::zero() /*controllare valore di verita in C */ {k=-k;}
                b=a;
                a=r;} else {a=a-b;}
    }
    }
pub fn even_part (mut a:Self)-> (Self,Self){
    /*even_part(a).0 is the even part of a
      even_part(a).1 is the grates exponent exp of 2 such that 2^exp divides a */
    
    let mut exp: Self=Self::zero();
    /*calulate the even part of b */
    while a%2==Self::zero() {
        exp+=Self::one();
        a=a/Self::from(2);                
    }
    (Self::pow(Self::from(2),exp),exp)
}
pub fn odd_part (a:Self)-> Self {
    let e:Self=Self::even_part(a).0;
    a/e
}
pub fn check_sqrt_mod_prime(a:Self,p:Self)->MathResult{
let e:Self=Self::even_part(p-1).1;
let q:Self=Self::odd_part(p-1);
if a%p==Self::zero() {return Some(Self::zero());}
/*find generator */
let mut rng = rand::thread_rng();
let mut n:u64=rng.gen();
n = loop {
    if Self::kronecker( n.into(), p)!=-1{n=rng.gen()} else {break n;}
};
let z:Self=Self::pow_mod(&n.into(), q,p);
/*initialize */
let mut y:Self=z.clone();
let mut r:Self=e.clone();
let mut x:Self=Self::pow_mod(&a, (q-1)/Self::from(2),p);
let mut b:Self=(((Self::pow_mod(&x, Self::from(2),p)))*a)%p;
x=(a*x)%p;

/*find exponent */
loop{
      let mut m:Self=Self::one();
            if b%p==Self::one() {break Some(x);} else {
              
                m=loop {
                        if Self::pow_mod(&b,Self::pow(Self::from(2), m),p)==Self::one()
                        {break m;} else{m+=Self::one();}
                        }};
            if m==r {break None;}

            let t:Self=Self::pow_mod(&y, Self::pow(Self::from(2), r-m-Self::one()),p);
            y=Self::pow_mod(&t, Self::from(2),p);
            r=m%p;
            x=(x*t)%p;
            b=(b*y)%p;
            
            }

}
pub fn pow_mod(&self,n: Self,p:Self) -> Self {
    let mut y:Self=Self::one();
    if n == Self::zero() {
        return y;}
    
        let mut aux:Self=n;
        let mut z=*self;
        
        loop{
        if aux%2==Self::one() {y=(z*y)%p}
      
        aux=aux/2;
       
        
        if aux==Self::zero() {break y;} else {z=(z*z)%p;}
        
        }
}
pub fn cubic_root_one_mod_p(p:Self)->RootsResult{
/*check if p is prime... */

/*if p is prime then try to solve z^3=1. z^3=1 is the trivial solution.
 z^3-1=0 iif (z-1)(z^2+z+1)=0. So non trivial solutions exist if delta=-3 is a quadratic residue.
 We chech if -3 is a quadratic residue witch check_sqrt_mod_prime. If true we solve the quadratic equation*/
let res:MathResult=Self::check_sqrt_mod_prime(p-3, p);
let delta:Self;
let z_1:Self;
let z_2:Self;
if res.is_some() {
    delta=res.unwrap();
    z_1= (((p-1)+delta)%p)*Self::pow_mod(&Self::from(2), p-2, p)%p;
    z_2=(((p-1)-delta)%p)*Self::pow_mod(&Self::from(2), p-2, p)%p;
        /*we put the smallest one for first and the biggest than */

    if z_1<z_2 {return [Some(z_1),Some(z_2)];} else {return [Some(z_2),Some(z_1)]}}
    else {return [None,None];};
    /*we put the smallest one for first and the biggest than */
   


}
pub fn expansion_base_c_root(mut n: Self,p:Self){
    /*This function takes as inputs n mod p 
 and if there exsists non trivial cubic roots, it outputs
 [t0,t1,t2] such that n=t0+t1z+t2z^2, where z is a non trivial cubic roots  */
let t0:Self;
let mut t1:Self=Self::zero();
let mut t2:Self=Self::zero();
let res:RootsResult=Self::cubic_root_one_mod_p(p);
let z:Self;
if res[0].is_some() {
    z=res[0].unwrap();
    let z_quad:Self=Self::pow_mod(&z,Self::from(2),p);
    while n>=z_quad{
        n=n-z_quad;
        t2+=Self::one();
    }
    while n>=z{
        n=n-z;
        t1+=Self::one();
    }
    t0=n;
    println!("Expantion in base cubic root is {:?}",[t0,t1,t2]);
    }
else {println!("There's no cubic root mod {}",p);};



}

}




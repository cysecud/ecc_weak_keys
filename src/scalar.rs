use uint::construct_uint;
use rand::Rng;
#[derive(Debug)]
pub enum MathError {
    QuadraticNonResidueModP
}
pub type MathResult = Option<U256>;
construct_uint!{
pub struct U256(8);
}
impl U256{
    pub fn random()->Self {
        let mut rng = rand::thread_rng();
        let mut v_rand = [0u64; 8];
        rng.fill(&mut v_rand);
        Self(v_rand)
    } 
    pub fn to_binary(mut n:Self)->Vec<Self>{
        let mut aux_bin: Vec<Self>=Vec::new();
        let bin:Vec<Self>;
        for _i in 0..n.bits(){
            aux_bin.push(n%Self::from(2));
            n=n/2;};
        bin= aux_bin.into_iter().rev().collect();
       bin}
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
        
        }/*now b must me odd */

    /*Recuce size once */
        a=a%(b);

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
    
    if a==Self::zero() {return Some(Self::zero());}
/*find generator */
    let mut rng = rand::thread_rng();
    let mut nrand:u8=rng.gen();
    let mut n:Self=Self::from(nrand)%p;
n = loop {
    if Self::kronecker( n, p)!=-1{
        nrand=rng.gen();
        n=Self::from(nrand)%p;} else {break n;}
};

    let z=Self::pow_mod(&n,q,p);
/*initialize */
    let mut y:Self=z.clone();
    let mut r:Self=e.clone();
    let mut x:Self=Self::pow_mod(&a, (q-Self::one())/Self::from(2),p);
    let mut b:Self=(Self::pow_mod(&mut x, Self::from(2),p)*a)%p;
    x=(a*x)%p;

/*find exponent */
loop{
    let mut m:Self=Self::one();
            if b==Self::one() {break Some(x);} else {
            
               m=loop { 
                if Self::pow_mod(& mut b, Self::pow(Self::from(2),m),p)==Self::one()
                        {break m;} else{m+=Self::one();}
                        }};
            if m==r {break None;}
            
            let mut t:Self=Self::pow_mod(&mut y, Self::pow(Self::from(2), r-m-Self::one()),p);
            
            y=Self::pow_mod(&mut t, Self::from(2),p);
            r=m%p;
            x=(x*t)%p;
            b=(b*y)%p;
            }

}

pub fn pow_mod(g:&Self,n: Self,p:Self) -> Self {
    let mut y:Self=Self::one();
    if n == Self::zero() {
        return y;}
    
        let mut aux:Self=n;
        let mut z:Self=g.clone();
        
        loop{
        if aux%2==Self::one() {y=(z*y)%p}
      
        aux=aux/2;
       
        
        if aux==Self::zero() {break y;} else {z=(z*z)%p;}
        
        }
}
}

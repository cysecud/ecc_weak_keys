use core::ops::{Add, Mul, Sub};
use crate::scalar::U256;
use std::collections::HashMap;
use reikna::totient::totient;
use reikna::factor::quick_factorize;
#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq,PartialOrd, Ord)]
pub struct FieldPoint {
    pub num: U256,
    pub prime: U256,
}
impl FieldPoint {
    pub fn new(num: U256, prime: U256) -> FieldPoint {
        if num > prime {
            panic!("Not a valid input for a field point, num should be nonnegative and less than prime, obtained {}", num);
        } else {
            FieldPoint {num:num, prime:prime}
        }
    }
    pub fn negate(self)->FieldPoint{
    FieldPoint {num:self.prime-self.num,
    		prime: self.prime}
    }
    pub fn multiple(&self,n:U256)->FieldPoint {
    FieldPoint{num:self.num*n%self.prime,
    		prime:self.prime}
    }
    pub fn inverse (& mut self)-> FieldPoint {
            FieldPoint::power(self, self.prime-U256::from(2))
    }
    pub fn inverse_mersenne(&self)->FieldPoint{
        /*Algorithm 11.9 Prime field inversion

        It's reported to be faster than the extended
        Euclidean algorithm for some types of primes, for example Mersenne primes.*/
        let mut z=*self;
        let mut u=FieldPoint::new(U256::one(),self.prime);
        let mut q:U256;
        while z.num!=U256::one() {
            q=self.prime/z.num;
            let num=(z.prime-q*z.num)%self.prime;
            z=FieldPoint{num:num,prime:self.prime};
            u=FieldPoint::multiple(&u, q).negate();
            //println!("z is {:?}, u is {:?}",z,u);
        }
        u
    }
    pub fn rand_mod(p:U256)-> FieldPoint{
    	let num=U256::random()%p;
    	FieldPoint::new(num,p)
    }

    pub fn square(&mut self)->Self {
        let num=U256::pow_mod(&self.num, U256::from(2), self.prime);
        FieldPoint::new(num, self.prime)
    }
    pub fn double(self)->Self{self+self}
    pub fn power(&mut self,n: U256) -> Self {
        /*Montgomery ladder */
        
        let bin=U256::to_binary(n);
        let mut p0=FieldPoint::new(U256::one(),self.prime);
        let mut p1=*self;
        for i in bin.iter() {
            if *i==U256::zero() {
                p1=p1*p0;
                p0=p0.square();
                
                /* println!("p0 is {:?}",p0);
                println!("p1 is {:?}",p1); */
            }
            else {
                p0=p0*p1;
                p1=p1.square();
              /*   println!("p0 is {:?}",p0);
                println!("p1 is {:?}",p1); */
                }
        }
        p0
        
    }
  }
/*-------------------Find a generator of Fp*--------------------*/
impl FieldPoint {
pub fn znprimroot(prime:&U256) ->  FieldPoint {

let  factor = quick_factorize((prime.as_u128()-1).try_into().unwrap());
//prendo i fattori primi senza molteplicità
	let mut freq_vec: HashMap<u64, u64> = HashMap::new();
	for i in &factor {
	let freq: &mut u64 = freq_vec.entry(*i).or_insert(0);
      	*freq += 1;
	};
let prime_factor = Vec::from_iter(freq_vec.keys());
//inizializazione
let mut alpha = U256::one();
/*
let mut z=FieldPoint::new(alpha,*prime);//alla fine dell'algoritmo questo sarà il generatore
*/

'outer:loop{
alpha+=U256::one();
let mut z=FieldPoint::new(alpha,*prime);
let mut i:usize=0;
	loop{
	let  zn=FieldPoint::power(&mut z,(*prime-1)/(*prime_factor[i]));

	if zn.num==U256::one() {break;} else {i+=1;};
	if i>=(prime_factor.len()).try_into().unwrap() {break 'outer z;};	
	}

     }
   }
}

impl FieldPoint {
pub fn znorder(g:&mut FieldPoint) -> U256 {
/*--------------------order of the integermod x in (Z/nZ)*----------------------*/	

/* Questa parte di algoritmo è generica; calcola l'ordine di
un elemento g in G data la cardinalità di G. 
Qui ci stiamo concentrando sul gruppo moltiplicativo Z/nZ*,
che ha cardinalità phi(n). Abbiamo dunque calcolato preliminarmente
tale valore phi(n) ed implementato l'algoritmo mediante la fattorizzazione
di phi(n).
Nel caso di un gruppo generico G useremo e=n=|G| */
	let n=U256::from(totient(g.prime.as_u64()));
	let factor = quick_factorize(n.as_u64());
	
	let mut freq_vec: HashMap<U256, U256> = HashMap::new();
	for i in &factor {		
		let freq: &mut U256 = freq_vec.entry((*i).try_into().unwrap()).or_insert(U256::zero());
		*freq += U256::one();
   	};
//vettore dei fattori primi di n senza ripetizioni 
let primes = Vec::from_iter(freq_vec.keys());

//vettore delle molteplicità dei fattori primi di n
let mut mol = Vec::new();
	for i in 0.. primes.len(){
	mol.push(factor.iter().filter(|&n| *n ==(primes[i].as_u64())).count())};
	
let mut e=n;

//let mut r=g.clone();
let mut i:usize = 0;

	loop{
		i+=1;
		if i>primes.len() {break e;}
		else {e=e/(primes[i-1].pow(U256::from(mol[i-1])));
		};
		let mut r=FieldPoint::power(g,e);
		while r.num != U256::one() {r=FieldPoint::power(&mut r,*primes[i-1]);
				  e=e*(*primes[i-1]);
				}
	}
}
/*-----------------Baby-Step-Giant-Step in Fp*-----------------------*/
pub fn bsgs(g:& mut FieldPoint,h:& mut FieldPoint,p:& U256) {
	let mut x=g.clone();
	let ord:U256= FieldPoint::znorder(&mut x);
	let  m =U256::integer_sqrt(&ord);

	let mut bs:HashMap<FieldPoint,U256>=HashMap::new();
	//these are the babysteps g^(i+1)
	for i in 0..m.as_u128() {
    		bs.insert(FieldPoint::power(g,U256::from(i+1)),U256::from(i));};
    	//calculate the inverse of the generator 
    	println!("bs list:{:?}",bs);	
	let mut inv_g=FieldPoint::power(g,p-2);
	let mut j:U256 = U256::zero();
	'outer: loop {
		for (k,v) in bs.iter() {
			//check if the babysteps in the hash match the giant-step g^(-m(j+1)), which is inv_g^(m(j+1))			
			
			if *k==	*h*FieldPoint::power(&mut inv_g,(j+1)*m)
				{let alfa=((*v+1)+m*(j+1))%p;
				println!("Found a match: alfa is {}",alfa); break 'outer;}}
			j=j+1;
			if j==m+1 {println!("no match found");
						break 'outer;}
			
		}
			                              
}
}
impl Add for FieldPoint{
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if self.prime == other.prime {
            FieldPoint {num: (self.num + other.num)%(self.prime), prime: self.prime}
        } else {
            panic!("Cannot add these field points, different prime values {},{}",self.prime,other.prime);
        }
    }}
impl Sub for FieldPoint {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        if self.prime == other.prime {
            FieldPoint {num: (self.num + other.negate().num)%(self.prime), prime: self.prime}
        } else {
            panic!("Cannot subtract these field points, different prime values {},{}",self.prime,other.prime);
        }
    }
}
impl Mul for FieldPoint {
    type Output = Self;
    #[inline]
    fn mul(self, other: Self) -> Self {
        if self.prime == other.prime {
            FieldPoint {num: (self.num*other.num)%(self.prime), prime: self.prime}
        } else {
            panic!("Cannot multiply these field points, different prime values, {},{}",self.prime,other.prime);
        }
    }
}


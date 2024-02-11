p=115792089237316195423570985008687907852837564279074904382605163141518161494337;

B=2^32;
/*Compute the divisors below the bound B*/
div=divisors(p-1);
div32=div;
foreach(div32,d,if(d>B,{div32=setminus(div32,[d])}));
s=0;
foreach(div32,d,s=s+eulerphi(d));
n32=log(s)/log(2); /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
Rp32=div32;
foreach(Rp32,d,Rp32=setminus(Rp32,setminus(divisors(d),[d])));
t=0;
foreach(Rp32,d,t=t+2*ceil(sqrt(d)));
c32=log(t)/log(2); /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

B=2^64;
/*Compute the divisors below the bound B*/
div=divisors(p-1);
div64=div;
foreach(div64,d,if(d>B,{div64=setminus(div64,[d])}));
s=0;
foreach(div64,d,s=s+eulerphi(d));
n64=log(s)/log(2); /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
Rp64=div64;
foreach(Rp64,d,Rp64=setminus(Rp64,setminus(divisors(d),[d])));
t=0;
foreach(Rp64,d,t=t+2*ceil(sqrt(d)));
c32=log(t)/log(2); /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

B=2^128;
/*Compute the divisors below the bound B*/
div=divisors(p-1);
div128=div;
foreach(div128,d,if(d>B,{div128=setminus(div128,[d])}));
s=0;
foreach(div128,d,s=s+eulerphi(d));
n128=log(s)/log(2); /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
Rp128=div128;
foreach(Rp128,d,Rp128=setminus(Rp128,setminus(divisors(d),[d])));
t=0;
foreach(Rp128,d,t=t+2*ceil(sqrt(d)));
c128=log(t)/log(2); /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

B=2^160;
/*Compute the divisors below the bound B*/
div=divisors(p-1);
div160=div;
foreach(div160,d,if(d>B,{div160=setminus(div160,[d])}));
s=0;
foreach(div160,d,s=s+eulerphi(d));
n160=log(s)/log(2); /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
Rp160=div160;
foreach(Rp160,d,Rp160=setminus(Rp160,setminus(divisors(d),[d])));
t=0;
foreach(Rp160,d,t=t+2*ceil(sqrt(d)));
c160=log(t)/log(2); /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

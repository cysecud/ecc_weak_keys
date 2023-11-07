/* Analysis of weak keys on the elliptic curve P-224 */
p=26959946667150639794667015087019625940457807714424391721682722368061;
B=2^120;
/*Compute the divisors below the bound B*/
div=divisors(p-1);
divB=div;
foreach(divB,d,if(d>B,{divB=setminus(divB,[d])}));
s=0;
foreach(divB,d,s=s+eulerphi(d));
nB=log(s)/log(2); /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
RpB=divB;
foreach(RpB,d,RpB=setminus(RpB,setminus(divisors(d),[d])));
t=0;
foreach(RpB,d,t=t+2*ceil(sqrt(d)));
cB=log(t)/log(2); /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

default(parisize, 120000000);
p=0xaadd9db8dbe9c48b3fd4e6ae33c9fc07cb308db3b3c9d20ed6639cca70330870553e5c414ca92619418661197fac10471db1d381085ddaddb58796829ca90069;

log(p)/log(2)
B=2^32;
div=divisors(p-1);
divB=div;
foreach(divB,d,if(d>B,{divB=setminus(divB,[d])}));
s=0;
foreach(divB,d,s=s+eulerphi(d));
print1("nB curve brainpoolP512r1, B=2^32")
nB=log(s)/log(2) /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
RpB=divB;
foreach(RpB,d,RpB=setminus(RpB,setminus(divisors(d),[d])));
t=0;
foreach(RpB,d,t=t+2*ceil(sqrt(d)));
print("cB curve brainpoolP512r1, B=2^32")
cB=log(t)/log(2) /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

B=2^64;

divB=div;
foreach(divB,d,if(d>B,{divB=setminus(divB,[d])}));
s=0;
foreach(divB,d,s=s+eulerphi(d));
print1("nB curve brainpoolP512r1, B=2^64")
nB=log(s)/log(2) /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
RpB=divB;
foreach(RpB,d,RpB=setminus(RpB,setminus(divisors(d),[d])));
t=0;
foreach(RpB,d,t=t+2*ceil(sqrt(d)));
print("cB curve brainpoolP512r1, B=2^64")
cB=log(t)/log(2) /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

B=2^128;

divB=div;
foreach(divB,d,if(d>B,{divB=setminus(divB,[d])}));
s=0;
foreach(divB,d,s=s+eulerphi(d));
print1("nB curve brainpoolP512r1, B=2^128")
nB=log(s)/log(2) /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
RpB=divB;
foreach(RpB,d,RpB=setminus(RpB,setminus(divisors(d),[d])));
t=0;
foreach(RpB,d,t=t+2*ceil(sqrt(d)));
print("cB curve brainpoolP512r1, B=2^128")
cB=log(t)/log(2) /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

B=2^160;

divB=div;
foreach(divB,d,if(d>B,{divB=setminus(divB,[d])}));
s=0;
foreach(divB,d,s=s+eulerphi(d));
print1("nB curve brainpoolP512r1, B=2^160")
nB=log(s)/log(2) /*this is the base-2 log of the number of weak keys with order bounded by B*/
/*We compute RpB, i.e., the set of divisors d1<d2<...<dk<=B, such that di is NOT a divisor of dj, for all 1<=i<j<=k */
RpB=divB;
foreach(RpB,d,RpB=setminus(RpB,setminus(divisors(d),[d])));
t=0;
foreach(RpB,d,t=t+2*ceil(sqrt(d)));
print("cB curve brainpoolP512r1, B=2^160")
cB=log(t)/log(2) /*this is the base-2 log of the worst-case number of elliptic-curve-group-operation 
required to test whether a key comes from a subgroup of order bounded by B*/

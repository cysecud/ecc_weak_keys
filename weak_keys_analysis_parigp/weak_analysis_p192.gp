p=6277101735386680763835789423176059013767194773182842284081;
B=2^32;
div=divisors(p-1);
divB=div;
foreach(divB,d,if(d>B,{divB=setminus(divB,[d])}));

/*divB=[1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 149, 192, 298, 447, 596, 631, 894, 1192, 1262, 1788, 1893, 2384, 2524, 3576, 3786, 4768, 5048, 7152, 7572, 9536, 10096, 14304, 15144, 20192, 28608, 30288, 40384, 60576, 94019, 121152, 188038, 282057, 376076, 564114, 752152, 1128228, 1504304, 2256456, 3008608, 4512912, 6017216, 9025824];*/ /*these are the divisors below the bound B*/
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

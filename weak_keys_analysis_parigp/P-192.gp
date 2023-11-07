q = 0xfffffffffffffffffffffffffffffffeffffffffffffffff
a = Mod(0xfffffffffffffffffffffffffffffffefffffffffffffffc, q)
b = Mod(0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1, q)
E = ellinit([a, b])
E[16][1] = 0xffffffffffffffffffffffff99def836146bc9b1b4d22831 * 0x1
G = [Mod(0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012, q), Mod(0x07192b95ffc8da78631011ed6b24cdd573f977a11e794811, q)]
p=ellorder(E,G)

alfa=3774452713541362385671898297224861821975924809030444427889;
p=ellorder(E,G);
Q=ellmul(E,G,alfa);
z=znprimroot(p);
bsgs(E,G,Q,d)=
{
m=ceil(sqrt(d))+1;
zd=lift(Mod(z,p)^((p-1)/d));
bs=Vec(0,m);
for(i=1,m,
	bs[i]=ellmul(E,Q,lift(Mod(zd,p)^i));
	for(j=1,i,
		if(bs[j]==ellmul(E,G,lift(Mod(zd,p)^(m*i))),
			
			print(lift(Mod(zd,p)^((m*i-j)%d)))
			) 
			)
	)
};

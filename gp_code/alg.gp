
print("Type key_pairs_gen() to generate your private key and public key");
print("Type test_key(you_public_key,bound_for_the_test) to test your key");
/* 
q = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f;
a = Mod(0x0000000000000000000000000000000000000000000000000000000000000000, q);
b = Mod(0x0000000000000000000000000000000000000000000000000000000000000007, q);
E = ellinit([a, b],q);
P = [Mod(55066263022277343669578718895168534326250603453777594175500187360389116729240,q),Mod(32670510020758816978083085130507043184471273380659243275938904335757337482424,q)];
p=115792089237316195423570985008687907852837564279074904382605163141518161494337;
z=Mod(7, p);
div32 = [18051648]; */
bsgs(public_key,d)=
{
	match=0; /*boolean value false if no match is found*/
	print("Point to be tested:\n " public_key);
	print("Divisor: " d);
m=ceil(sqrt(d))+1;
zd=lift(Mod(z,p)^((p-1)/d));
bs=Vec(0,m);
for(i=1,m,
	bs[i]=ellmul(E,public_key,lift(Mod(zd,p)^i)));

for(i=1,m,
		gian_step=ellmul(E,P,lift(Mod(zd,p)^(m*i)));
		foreach(bs,baby_step,
			if (baby_step==gian_step,
				print("Match found!");
				j=select((x) -> x == baby_step, bs, 1);
				private_key=lift(Mod(zd,p)^((m*i-j[1])%d));
				print("private key detected:\n" private_key);
				match=1;
				break;
				)
				)
			
	);
if (match == 0,print("No match found! Private key is safe within this bound!"));

};
test_key(public_key,bound)={
	
	if(bound==32, 
		foreach(div32,d,bsgs(public_key,d)),
		if(bound==64, 
			foreach(div64,d,bsgs(public_key,d)),
				if(bound==128, 
					foreach(div128,d,bsgs(public_key,d)),
					if(bound==160, 
					foreach(div160,d,bsgs(public_key,d)))
				)
		), print("Wrong bound size. Available size are 32,64,128,160."))
};
key_pairs_gen()={
	d=div32[1];
	zd=115481771728459905245102424859900657047113141323743738905491223467302634970004;
	sk=lift(Mod(zd,p)^random(d));
	pk=ellmul(E,P,sk);
	print("Your private and public keys have been generated!\nType sk to visualize you secret key and pk to visualize the public key.");

}

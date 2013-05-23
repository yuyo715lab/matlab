function YprimeOpen = open(from,to,P,Q,RHO,GorL)
	[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,GorL] = Define();
	R(from,to) = 2 * R(from,to);
	R(to,from) = R(from,to);
	Y = admittance(N,R,Tr);
	[YprimeOpen,numG,numL] = Yprime(N,Y,Ps,Qs,PQorPV,P,Q,RHO,GorL);
end

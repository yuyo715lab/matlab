function YprimeOpen = open(from,to,P,Q,RHO,GorL)
	[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,GorL] = Define();
	R(from,to) = 2 * R(from,to);
	R(to,from) = R(from,to);
	R(5,5) = 1/(0.088i + 0.0765i);
	R(7,7) = 1/(0.0765i + 0.0745i);
	Y = admittance(N,R,Tr);
	[YprimeOpen,numG,numL] = Yprime(N,Y,Ps,Qs,PQorPV,P,Q,RHO,GorL);
end

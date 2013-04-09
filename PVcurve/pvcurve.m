%============= Define.m ================
[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ] = Define();
%============= Define.m ================

%============= admittance.m ================
Y = admittance(N,R,Tr);
%============= admittance.m ================


[P,Q,RHO,THEATA] = main(Y);
Y0 = Y;


%######### specified node #########
SN = 5;
%######### specified node #########

%######### Cycle Number #########
CN = 20;
%######### Cycle Number #########

RHO(SN)
-P(SN)

for k = 1:CN
  Y = Y0;
  Y(SN,SN) = Y0(SN,SN) + (0.8 + 0.3*(k-1))*(-P(SN) - Q(SN)*i) / (RHO(SN)^2);
  [P,Q,RHO,THEATA] = ...
      main(Y);
  RHO(SN) 
 - P(SN)
end

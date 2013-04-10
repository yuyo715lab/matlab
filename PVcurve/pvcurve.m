%######### specified node #########
SN = 5;
%######### specified node #########


%######### First PQ-proportion #########
fprop = 0;
%######### First PQ-proportion #########

%######### Delta propotion #########
dprop = 0.01;
%######### Delta proportion #########

%######### P limit #########
Plimit = 4.06;
%######### P limit #########

%============= Define.m ================
[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,Pss] = Define(1,SN);
%============= Define.m ================

CN = 0;
prop = fprop;
while -Pss * prop < Plimit
  CN = CN + 1;
  prop = prop + dprop;
end


XY = zeros(2,CN);

k = 1;
prop = fprop;
while -Pss * prop < Plimit

  %============= Define.m ================
  [N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,Pss] = Define(prop,SN);
  %============= Define.m ================

  [P,Q,RHO,THEATA] = ...
      main(N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ);
  XY(1,k) = -P(SN);
  XY(2,k) = RHO(SN);
  prop = prop + dprop;
  k = k + 1;
end

plot(XY(1,:),XY(2,:))








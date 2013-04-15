clear all

startT = clock();
startCpuT = cputime;

%######### specified node #########
SN = 8;
%######### specified node #########

%######### initial V #########
ei = 1.0;
%######### initial V #########

%######### First PQ-proportion #########
fprop = 0;
%######### First PQ-proportion #########

%######### Delta propotion #########
dprop = 0.0001;
%######### Delta proportion #########

%######### P limit #########
Plimit = 4.6724;
%######### P limit #########

%============= Define.m ================
[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,Pss] = Define(1,SN,ei);
%============= Define.m ================

CN = 0;
prop = fprop;
while -Pss * prop < Plimit
  CN = CN + 1;
  prop = prop + dprop;
end


XY = zeros(2*CN,2);

k = 1;
prop = fprop;
while -Pss * prop < Plimit

  %============= Define.m ================
  [N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,Pss] = Define(prop,SN,ei);
  %============= Define.m ================

  [P,Q,RHO,THEATA] = ...
      main(N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ);
  XY(k,1) = -P(SN);
  XY(k,2) = RHO(SN);
  prop = prop + dprop;
  k = k + 1;
end

ei = 0;
prop = fprop;
while -Pss * prop < Plimit

  %============= Define.m ================
  [N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ,Pss] = Define(prop,SN,ei);
  %============= Define.m ================

  [P,Q,RHO,THEATA] = ...
      main(N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ);
  XY(k,1) = -P(SN);
  XY(k,2) = RHO(SN);
  prop = prop + dprop;
  k = k + 1;
end





ntime=cputime-startCpuT;
nhour = floor(ntime/60/60);
nmin = floor((ntime-nhour*3600)/60);
nsec = ntime-nhour*3600-nmin*60;
disp(sprintf('%s%s', 'Start time',datestr(startT,31)));
disp(sprintf('%s%s', 'End time',datestr(clock,31)));
disp(sprintf('%s%d%s%02d%s%04.1f%s', 'time required',nhour,'h',nmin,'m',nsec,'s'));

%plot(XY(1,:),XY(2,:))
scatter(XY(:,1),XY(:,2))







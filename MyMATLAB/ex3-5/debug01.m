clear all


%============= Define.m ================
[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ] = Define();
%============= Define.m ================

%============= admittance.m ================
Y = admittance(N,R,Tr);
%============= admittance.m ================

%============= debug_PandQ.m ================
[dP,dQ,dV] = debug_PandQ(N,Y,PQ,Ps,Qs,Vs,e,f);
%============= debug_PandQ.m ================

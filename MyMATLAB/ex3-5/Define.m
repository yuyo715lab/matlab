function [N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ] = Define()


%%%%%%%%%%% Bus Number %%%%%%%%%%%
N = 4;
%%%%%%%%%%% Bus Number %%%%%%%%%%%

%%%%%%%%% Reference node %%%%%%%%%
Ref = 1;
%%%%%%%%% Reference node %%%%%%%%%

%%%%%%%%% PQ or PV %%%%%%%%
%% Reference node - 0 %% PQ node - 1 %% PV node - 2 %% 
PQorPV = zeros(1,N);
PQorPV = [0 1 1 2];
NonRef = find(PQorPV);
%%%%%%%%% PQ or PV %%%%%%%%%


%%%%%%%%%%% Impedance %%%%%%%%%%%
R = zeros(N,N);

R(1,1) = -30i;
R(1,2) = 0.08 + 0.4i;
R(1,3) = 0.12 + 0.5i;
R(1,4) = 0;
R(2,1) = R(1,2);
R(2,2) = -34i;
R(2,3) = 0.1 + 0.4i;
R(2,4) = 0;
R(3,1) = R(1,3);
R(3,2) = R(2,3);
R(3,3) = -29i;
R(3,4) = 0.3i;
R(4,1) = R(1,4);
R(4,2) = R(2,4);
R(4,3) = R(3,4);
R(4,4) = 0;


%%%%%%%%%%% Impedance %%%%%%%%%%%


%%%%%%%%%%% Off norminal turn ratio %%%%%%%%%%%
Tr = zeros(N,N);

Tr(1,1) = 0;
Tr(1,2) = 0;
Tr(1,3) = 0;
Tr(1,4) = 0;
Tr(2,1) = 0;
Tr(2,2) = 0;
Tr(2,3) = 0;
Tr(2,4) = 0;
Tr(3,1) = 0;
Tr(3,2) = 0;
Tr(3,3) = 0;
Tr(3,4) = 1.1;
Tr(4,1) = 0;
Tr(4,2) = 0;
Tr(4,3) = 0;
Tr(4,4) = 0;

%%%%%%%%%%% Off norminal turn ratio %%%%%%%%%%%


%%%%%%%%% Node Voltage %%%%%%%%%
e = zeros(1,N);
f = zeros(1,N);
Vs = zeros(1,N);
V = zeros(1,N);
dV = zeros(1,N);

e = [1.05 1.0 1.0 1.0];
Vs(4) = 1.1;

%%%%%%%%% Node Voltage %%%%%%%%%


%%%%%%%%% P and Q %%%%%%%%%
Ps = zeros(1,N);
Qs = zeros(1,N);
PQ = zeros(1,N);

Ps = [0 -0.55 -0.3 0.5];
Qs = [0 -0.13 -0.18 0];
%%%%%%%%% P and Q %%%%%%%%%




clear all




%--------- power flow calculation ---------
[N,Ps,Qs,PQorPV,P,Q,RHO,THEATA,Y,GorL] = PowerCalc();
%--------- power flow calculation ---------

%--------- Yprime ---------
[Yprime,numG,numL] = Yprime(N,Y,Ps,Qs,PQorPV,P,Q,RHO,GorL);
%--------- Yprime ---------



%--------- YprimeEF ---------
EarthFault = 5;
[YprimeEF] = YprimeEF(N,Y,Ps,Qs,PQorPV,P,Q,RHO,GorL,EarthFault);
%--------- YprimeEF ---------

Glabel = zeros(1,numG);
Glabel = find(GorL == 0);

%--------- eqipment constant ---------
[xd,xdd,xddd,xq,xqq,xqqq,xl,Td,Tdd,Tq,Tqq,R,Kg,Tg,Ka,Ta,D,H] ...
      = equipment(numG);
%--------- eqipment constant ---------

%--------- state ---------
[Idq,Edq] = state(N,numG,GorL,RHO,THEATA,Yprime,Glabel,R,xq);
%--------- state ---------
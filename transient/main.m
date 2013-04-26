
clear all


max = 200;

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
[xd,xdd,xddd,xq,xqq,xqqq,xl,Td,Tdd,Tq,Tqq,Rg,KG,TG,KA,TA,D,H,Kd,Kq] ...
      = equipment(numG);
%--------- eqipment constant ---------

%--------- state ---------
[Idq,Vdq,Edq,deltaEq,id,iq,vd,vq,ef0,Pe,eq,eqq,ed,edd,egd,egq,Egd,Egq] = ...
    state(N,numG,GorL,RHO,THEATA,Yprime,Glabel,xd,xq,Rg,Kd,Kq,xdd,xl,xqq,xqqq,xddd,max);
%--------- state ---------

%--------- Yg ---------
[Yg] = YG(numG,xddd,xqqq,Rg,deltaEq,YprimeEF);
%--------- Yg ---------

Egdq = zeros(numG*2,1);
for k = 1:numG
  Egdq(2*k-1) = Egd(k);
  Egdq(2*k) = Egq(k);
end

Idq = Yg * Egdq;
Vdq = inv(YprimeEF) * Idq;

%--------- id iq ---------
for k = 1:numG
  id(k) = real((Idq(2*k-1)+Idq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
  iq(k) = imag((Idq(2*k-1)+Idq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
  vd(k) = real((Vdq(2*k-1)+Vdq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
  vq(k) = imag((Vdq(2*k-1)+Vdq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
end
%--------- id iq ---------

%--------- Pe ---------
Pe(1,:) = vd .* id + vq .* iq + Rg .* (id.^2 + iq.^2);
%--------- Pe ---------


Runge_Kutta(P,numG,Pe,H,D,TG,KG,Td,xd,xdd,xddd,xl,id,Kd,Kd,Kq);
clear all

%for EarthFault = 1:9
max = 3000;


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
[Idq,Vdq,Edq,deltaEq,id,iq,vd0,vq0,ef0,Pe,eq,eqq,ed,edd,egd,egq,Egd,Egq] = ...
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

%eft_min = 17;
%eft_max = 45;
delta_for_plot = zeros(max,6);
w_for_plot = zeros(max,6);
v_for_plot = zeros(max,6);
for k = 40:50
	EarthFaultTime = k;
	[delta_for_plot,w_for_plot,v_for_plot] =...
		 Runge_Kutta2(P,numG,Pe,H,D,TG,KG,Td,Tdd,Tq,xd,xdd,xddd,xl,id,Kd,Kq,vd,vq,KA,TA,...
		 xq,xqq,xqqq,iq,Tqq,ef0,deltaEq,eq,eqq,ed,edd,vd0,vq0,Yg,YprimeEF,max,Yprime,...
		 Rg,Glabel,EarthFaultTime,delta_for_plot,w_for_plot,v_for_plot);
% $$$ 	if abs((delta(1,2) - delta(1,1)) - (delta(max,2)-delta(max,1)))/pi*180 < 0.1 &&...
% $$$ 			 abs((delta(1,3) - delta(1,1)) - (delta(max,3)-delta(max,1)))/pi*180 < 0.1
% $$$ %		check(EarthFaultTime-19,1) = EarthFaultTime;
% $$$ %		check(EarthFaultTime-19,EarthFault+1)= 0;
% $$$ 	else							
		%		check(EarthFaultTime-19,1)= EarthFaultTime;
		%		check(EarthFaultTime-19,EarthFault+1)= 1;
% $$$ 		if abs((delta(1,2) - delta(1,1)) - (delta(max,2)-delta(max,1)))/pi*180 > 0.1 ||...
% $$$ 				 abs((delta(1,3) - delta(1,1)) - (delta(max,3)-delta(max,1)))/pi*180 > 0.1
% $$$ 			EarthFault
% $$$ 			EarthFaultTime
% $$$ 			break
% $$$ 		end							
%	EarthFaultTime				
end
%check
%clear all
%end

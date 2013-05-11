function [delta_for_plot,w_for_plot,v_for_plot] =Runge_Kutta3(P,numG,Pe,H,D,TG,KG,Td,Tdd,Tq,xd,xdd,xddd,xl,id,Kd,Kq,vd,vq,KA,TA,...
    xq,xqq,xqqq,iq,Tqq,ef0,deltaEq,eq_state,eqq_state,ed_state,edd_state,vd0,vq0,Yg,YprimeEF,...
    max,Yprime,Rg,Glabel,EarthFaultTime,delta_for_plot,w_for_plot,v_for_plot,dt,now_step,endTime)

f = 50; % frequency
w0 = 2*pi*f; % angular velocity

%$$$$$$$$$ initialization $$$$$$$$$
w = zeros(max,numG);
Pm = zeros(max,numG);
Pm0 = zeros(1,numG);
eq = zeros(max,numG);
ef = zeros(max,numG);
eqq = zeros(max,numG);
%eqqq = zeros(max,numG);
ed = zeros(max,numG);
edd = zeros(max,numG);
%eddd = zeros(max,numG);
V = zeros(max,numG);
delta = zeros(max,numG);
egd = zeros(max,numG);
egq = zeros(max,numG);
%$$$$$$$$$ initialization $$$$$$$$$

%((((((((( initial value )))))))))
delta(1,:) = deltaEq;
eq(1,:) = eq_state;
eqq(1,:) = eqq_state;
ed(1,:) = ed_state;
edd(1,:) = edd_state;
ef(1,:) = ef0;
for k = 1:numG
  w(1,k) = w0;
  Pm0(1,k) = P(Glabel(k));
  Pm(1,k) = P(Glabel(k));
  V(1,k) = sqrt(vd0(k)^2 + vq0(k)^2);
end
%((((((((( initial value )))))))))

%/////////////////////////////////////////////
%//////////// for loop ///////////////////////
for n = 1:max-1
  %<<<<<<<<<<<<<< earth fault >>>>>>>>>>>>>>>
  if n == EarthFaultTime/dt % earth fault occur at n*dt
    [Yg] = YG(numG,xddd,xqqq,Rg,deltaEq,Yprime);
    YprimeEF = Yprime;
  end
  %<<<<<<<<<<<<<< earth fault >>>>>>>>>>>>>>>
  delta_k = zeros(5,numG);
  w_k = zeros(5,numG);
  eq_k = zeros(5,numG);
  eqq_k = zeros(5,numG);
  ed_k = zeros(5,numG);
  edd_k = zeros(5,numG);
  ef_k = zeros(5,numG);
  Pm_k = zeros(5,numG);
  dt_k_tmp = [0;dt/2;dt/2;dt];
  dt_k = dt_k_tmp*ones(1,3);
  idqvdq_atk = zeros(20,numG);
%  id_atk = zeros(5,numG);
%  iq_atk = zeros(5,numG);
%  vd_atk = zeros(5,numG);
%  vq_atk = zeros(5,numG);
  for k = 1:4
    rk_egdatk = ...
	RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:));
    rk_egqatk = RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt_k(k,:).*eqq_k(k,:),eq(n,:)+dt_k(k,:).*eq_k(k,:));
    rk_egdqatk = RK_EGDQ(numG,delta(n,:)+dt_k(k,:).*delta_k(k,:),...
	rk_egdatk,rk_egqatk);
[id_atk iq_atk] = RK_idq(delta(n,:)+dt_k(k,:).*delta_k(k,:),Yg,...
	rk_egdqatk,numG);
[vd_atk vq_atk] = RK_vdq(delta(n,:)+dt_k(k,:).*delta_k(k,:),Yg,...
	rk_egdqatk,numG,YprimeEF);
%[vd_atk vq_atk] = RK_vdq(rk_egqatk,rk_egdatk,w0,xqqq,xddd,iq_atk,id_atk,Rg);
    idqvdq_atk(k+1,:) = id_atk;
	 idqvdq_atk(5+k+1,:) = iq_atk;
    idqvdq_atk(10+k+1,:) = vd_atk;
    idqvdq_atk(15+k+1,:) = vq_atk;
    delta_k(k+1,:) = RK_delta(w(n,:)+dt_k(k,:).*delta_k(k,:),dt,f,n);
    w_k(k+1,:) = RK_w(w(n,:)+dt_k(k,:).*w_k(k,:),dt,n,f,H,Pm(n,:)+dt_k(k,:).*Pm_k(k,:),...
	RK_Pe(rk_egdatk,idqvdq_atk(k+1,:),rk_egqatk,idqvdq_atk(5+k+1,:),w0,xddd,xqqq,w(n,:)+dt_k(k,:).*w_k(k,:)),D,numG);
    eq_k(k+1,:) = RK_eq(numG,Td,ef(n,:)+dt_k(k,:).*ef_k(k,:),xd,xdd,xddd,xl,Kd,...
	eqq(n,:)+dt_k(k,:).*eqq_k(k,:),...
	eq(n,:)+dt_k(k,:).*eq_k(k,:),w(n,:)+dt_k(k,:).*w_k(k,:),idqvdq_atk(k+1,:),n,w0);
    eqq_k(k+1,:) =RK_eqq(numG,Tdd,Kd,eq(n,:)+dt_k(k,:).*eq_k(k,:),eqq(n,:)+dt_k(k,:).*eqq_k(k,:),...
	w(n,:)+dt_k(k,:).*w_k(k,:),xdd,xl,idqvdq_atk(k+1,:),n,w0);
    ed_k(k+1,:) = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd(n,:)+dt_k(k,:).*edd_k(k,:),...
	w(n,:)+dt_k(k,:).*w_k(k,:),...
	idqvdq_atk(5+k+1,:),n,ed(n,:)+dt_k(k,:).*ed_k(k,:),w0);
    edd_k(k+1,:)=RK_edd(numG,Tqq,Kq,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:),...
	w(n,:)+dt_k(k,:).*w_k(k,:),xqq,xl,idqvdq_atk(5+k+1,:),n,w0);
    ef_k(k+1,:) = RK_ef(TA,KA,idqvdq_atk(10+k+1,:),idqvdq_atk(15+k+1,:),ef0,V,n,ef(n,:)+dt_k(k,:).*ef_k(k,:));
    Pm_k(k+1,:)=RK_Pm(Pm(n,:)+dt_k(k,:).*Pm_k(k,:),w(n,:)+dt_k(k,:).*w_k(k,:),numG,TG,KG,w0,n,Pm0);
  end
  
  delta(n+1,:) = delta(n,:) + dt./6.*(delta_k(2,:)+2.*delta_k(3,:)+2.*delta_k(4,:)+delta_k(5,:));
  w(n+1,:) = w(n,:) + dt./6.*(w_k(2,:)+2.*w_k(3,:)+2.*w_k(4,:)+w_k(5,:));
  eq(n+1,:) = eq(n,:) + dt./6.*(eq_k(2,:)+2.*eq_k(3,:)+2.*eq_k(4,:)+eq_k(5,:));
  eqq(n+1,:) = eqq(n,:) + dt./6.*(eqq_k(2,:)+2.*eqq_k(3,:)+2.*eqq_k(4,:)+eqq_k(5,:));
  ed(n+1,:) = ed(n,:) + dt./6.*(ed_k(2,:)+2.*ed_k(3,:)+2.*ed_k(4,:)+ed_k(5,:));
  edd(n+1,:) = edd(n,:) + dt./6.*(edd_k(2,:)+2.*edd_k(3,:)+2.*edd_k(4,:)+edd_k(5,:));
  ed(n+1,:) = ed(n,:) + dt./6.*(ed_k(2,:)+2.*ed_k(3,:)+2.*ed_k(4,:)+ed_k(5,:));
  ef(n+1,:) = ef(n,:) + dt./6.*(ef_k(2,:)+2.*ef_k(3,:)+2.*ef_k(4,:)+ef_k(5,:));
  Pm(n+1,:) = Pm(n,:) + dt./6.*(Pm_k(2,:)+2.*Pm_k(3,:)+2.*Pm_k(4,:)+Pm_k(5,:));

  egd(n+1,:) = RK_egd(Kq,xqq,xqqq,xl,edd(n+1,:),ed(n+1,:));
  egq(n+1,:) = ...
      RK_egq(Kd,xdd,xddd,xl,eqq(n+1,:),eq(n+1,:));
  rk_egdq = RK_EGDQ(numG,delta(n+1,:),egd(n+1,:),egq(n+1,:));
  [id_atn1 iq_atn1] = RK_idq(delta(n+1,:),Yg,rk_egdq,numG);
  [vd_atn1 vq_atn1] = RK_vdq(delta(n+1,:),Yg,rk_egdq,numG,YprimeEF);
%[vd_atn1 vq_atn1] = RK_vdq(egd(n+1,:),egq(n+1,:),w0,xqqq,xddd,iq_atn1,id_atn1,Rg);
  id(n+1,:) = id_atn1;
  iq(n+1,:) = iq_atn1;
  vd(n+1,:) = vd_atn1;
  vq(n+1,:) = vq_atn1;
  Pe(n+1,:) = RK_Pe(egd(n+1,:),id(n+1,:),egq(n+1,:),iq(n+1,:),w0,xddd,xqqq,w(n+1,:));


%  k = round(EarthFaultTime*100-37);
%  delta_for_plot(n,1) = dt*(n-1);	 
  %		w_for_plot(n,1) = dt*(n-1);
  %		v_for_plot(n,1) = dt*(n-1);
  
%  delta_for_plot(n,now_step) = (delta(n,2)-delta(n,1))/pi*180;
  %		delta_for_plot(n,2*(k-1)+1) = (delta(n,3) - delta(n,1))/pi*180;

  %		w_for_plot(n,k) = w(n,1)/2/pi;
  %		w_for_plot(n,3) = w(n,2)/2/pi;
  %		w_for_plot(n,4) = w(n,3)/2/pi;

  %v_for_plot(n,k) = sqrt(vd(n,1)^2+vq(n,1)^2);
  %		v_for_plot(n,3) = sqrt(vd(n,2)^2+vq(n,2)^2);
  %		v_for_plot(n,4) = sqrt(vd(n,3)^2+vq(n,3)^2);

  if mod(n,10000) == 0
    str = ['EFT ',num2str(EarthFaultTime),' n ',num2str(n)];
    disp(str)
  end
end
delta_for_plot(:,1) = [0:dt:endTime-dt]';
delta_for_plot(:,now_step) = (delta(:,2)-delta(:,1))/pi*180;
%//////////////// for loop ////////////////////
%/////////////////////////////////////////////
% $$$ 	if EarthFaultTime == 45
% $$$ 		plot(delta_for_plot(:,1),delta_for_plot(:,2))
% $$$ 		hold all
% $$$ 		plot(delta_for_plot(:,1),delta_for_plot(:,3))
% $$$ 	end
% $$$ 	hold all
% $$$ 	plot(v_for_plot(:,1),v_for_plot(:,3))
% $$$ 	hold all
% $$$ 	plot(v_for_plot(:,1),v_for_plot(:,4))

%	csvwrite('w_eftime.csv',w_for_plot);
%	csvwrite('v_eftime_001.csv',v_for_plot);
end


%////////////////////////// differential equation //////////////////////////
%--------- (8.2) ---------
function [ddeltadt] = RK_delta(w,dt,f,n)
w0 = 2*pi*f;
ddeltadt = w - w0;
end
%--------- (8.2) ---------

%--------- (8.3) ---------
function [dwdt] = RK_w(w,dt,n,f,H,Pm,Pe,D,numG)
dwdt = zeros(1,numG);
w0 = 2*pi*f;
dwdt = (w0./2./H.*(w0./w.*Pm - w0./w.*Pe - ...
    D.*(w./w0 - 1))) ;
end
%--------- (8.3) ---------

%--------- (8.4) ---------
function [deqdt] = RK_eq(numG,Td,ef,xd,xdd,xddd,xl,Kd,eqq,eq,w,id,n,w0)
deqdt = zeros(1,numG);
deqdt = 1./Td.*(ef + (xd - xdd).*(xdd - xddd)./((xdd - ...
    xl).^2).*Kd.*eqq - ...
    (1+(xd-xdd).*(xdd-xddd)./((xdd-xl).^2)).*eq - ...
    w./w0.*(xd-xdd).*(xddd-xl)./(xdd-xl).*id);
end
%--------- (8.4) ---------

%--------- (8.5) ---------
function [deqqdt] = RK_eqq(numG,Tdd,Kd,eq,eqq,w,xdd,xl,id,n,w0)
deqqdt = zeros(1,numG);
deqqdt = ...
    -1./Tdd./Kd.*(Kd.*eqq-eq+w./w0.*(xdd-xl).*id);
end
%--------- (8.5) ---------

%--------- (8.6) ---------
function [deddt] = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd,w,iq,n,ed,w0)
deddt = zeros(1,numG);
deddt = -1./Tq.*(-(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2).*Kq.*edd...
    +(1+(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2)).*ed ...
    -w./w0.*(xq-xqq).*(xqqq-xl)./(xqq-xl).*iq);
end
%--------- (8.6) ---------
%--------- (8.7) ---------
function [dedddt] = RK_edd(numG,Tqq,Kq,edd,ed,w,xqq,xl,iq,n,w0)
dedddt = zeros(1,numG);
dedddt = -1./Tqq./Kq.*(Kq.*edd-ed-w./w0.*(xqq-xl).*iq);
end
%--------- (8.7) ---------

%--------- (8.115) ---------
function [defdt] = RK_ef(TA,KA,vd,vq,ef0,V,n,ef)
defdt = -1./TA.*(ef+KA.*sqrt(vd.^2+vq.^2)...
    -(ef0+KA.*V(1,:))); %V hanyousei
end
%--------- (8.115) ---------

%--------- (8.116) ---------
function [dPmdt] = RK_Pm(Pm,w,numG,TG,KG,w0,n,Pm0)
dPmdt = zeros(1,numG);
dPmdt = -1./TG .* (Pm + KG ./ w0 .* w - (Pm0(1,:) + KG));
end
%--------- (8.116) ---------
%////////////////////////// differential equation //////////////////////////


function [id,iq] = RK_idq(delta,Yg,EGDQ,numG)
IGDQ = Yg*EGDQ;
id = zeros(1,numG);
iq = zeros(1,numG);
idq = zeros(numG,1);
idq = (IGDQ(1:2:end) + 1i.*IGDQ(2:2:end)) .* exp(1i.*(pi/2-delta.'));
id = real(idq.');
iq = imag(idq.');
end

function [vd,vq] = RK_vdq(delta,Yg,EGDQ,numG,YprimeEF)
VGDQ = YprimeEF\(Yg*EGDQ);
vd = zeros(1,numG);
vq = zeros(1,numG);
vdq = zeros(numG,1);
vdq = (VGDQ(1:2:end) + 1i.*VGDQ(2:2:end)) .* exp(1i*(pi/2-delta.'));
vd = real(vdq.');
vq = imag(vdq.');
end

function [EGDQ] = RK_EGDQ(numG,delta,egd,egq)
EGDQ = zeros(2*numG,1);
EGD = zeros(numG,1);
EGQ = zeros(numG,1);
EGD = egd .* sin(delta) + egq .* cos(delta);
EGQ = -egd .* cos(delta) + egq .* sin(delta);
EGDQ(1:2:end) = EGD;
EGDQ(2:2:end) = EGQ;
end

function [egd] = RK_egd(Kq,xqq,xqqq,xl,edd,ed)
egd = Kq.*(xqq-xqqq)./(xqq-xl).*edd + (xqqq-xl)./(xqq-xl).*ed;
end

function [egq] = RK_egq(Kd,xdd,xddd,xl,eqq,eq)
egq = Kd.*(xdd-xddd)./(xdd-xl).*eqq + (xddd-xl)./(xdd-xl).*eq;
end

function [Pe] = RK_Pe(egd,id,egq,iq,w0,xddd,xqqq,w)
Pe = (egd.*id+egq.*iq-(xddd-xqqq).*id.*iq).*w./w0;
end
%{
function [vd,vq] = RK_vdq(egd,egq,w0,xqqq,xddd,iq,id,R)
	vd = egd + w0 .* xqqq .*iq - R .* id;
	vq = egq + w0 .* xddd .*id - R .* iq;
end
%}

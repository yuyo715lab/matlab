function [] =Runge_Kutta(P,numG,Pe,H,D,TG,KG,Td,Tdd,Tq,xd,xdd,xddd,xl,id,Kd,Kq,vd,vq,KA,TA,...
		 xq,xqq,xqqq,iq,Tqq,ef0,deltaEq,eq_state,eqq_state,ed_state,edd_state);


%t(n+1) = t(n) + dt;
max = 200;
f = 50;

dt = 0.01;
w0 = 2*pi*f;
w = zeros(max,numG);
Pm = zeros(max,numG);
eq = zeros(max,numG);
ef = zeros(max,numG);
eqq = zeros(max,numG);
eqqq = zeros(max,numG);
ed = zeros(max,numG);
edd = zeros(max,numG);
eddd = zeros(max,numG);
V = zeros(max,numG);
delta = zeros(max,numG);

%((((((((( initial value )))))))))
delta(1,:) = deltaEq(1,:);
eq(1,:) = eq_state(1,:);
eqq(1,:) = eqq_state(1,:);
ed(1,:) = ed_state(1,:);
edd(1,:) = edd_state(1,:);
ef(1,:) = ef0;
for k = 1:numG
  w(1,k) = w0;
  Pm(1,k) = P(k);
  V(1,k) = sqrt(vd(k)^2 + vq(k)^2);
end
%((((((((( initial value )))))))))
n = 1;

%--------- k1 ---------
delta_k1 = RK_delta(w(n,:),dt,f,n)
w_k1 = RK_w(w(n,:),dt,n,f,H,Pm,Pe,D,numG)
eq_k1 = RK_eq(numG,Td,ef(n,:),xd,xdd,xddd,xl,...
    Kd,eqq,eq,w(n,:),id(n,:),n)
eqq_k1 = RK_eqq(numG,Tdd,Kd,eq(n,:),eqq(n,:),w(n,:),xdd,xl,id(n,:),n)
ed_k1 = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd,w(n,:),iq,n,ed(n,:))
edd_k1 = RK_edd(numG,Tqq,Kq,edd(n,:),ed(n,:),w(n,:),xqq,xl,iq,n)
ef_k1 = RK_ef(TA,KA,vd,vq,ef0,V,n,ef(n,:))
Pm_k1 = RK_Pm(Pm(n,:),w(n,:),numG,TG,KG,w0,n) 
%--------- k1 ---------
%--------- k2 ---------
wk_2 = RK_w(w(n,:) + dt .* w_k1./2 ,dt,n,f,H,Pm(n,:) + dt .* ...
    Pm_k1./2,Pe,D,numG)
Pm_k2 = RK_Pm(Pm(n,:)+dt.*Pm_k1./2,w(n,:)+dt.*w_k1./2,numG,TG,KG,w0,n)
%eq_k2 = RK_eq(numG,Td,ef(n,:)+ef_k1.*dt./2,xd,xdd,xddd,xl,Kd,...
 %   eqq(n,:)+eqq_k1.*dt./2,eq(n,:)+eq_k1.*dt./2,w(n,:)+w_k1.*dt./2,id(n,:)+id_k1.*dt./2,n)
%eqq_k2 = RK_eqq(numG,Tdd,Kd,eq(n,:)+eqq_k1.*dt./2,w(n,:)+w_k1.*dt./2,xdd,xl,id(n,:)+
%--------- k2 ---------
%k3 = RK_w(x(n) + dt/2 * k2);
%k4 = RK_w(x(n) + dt * k3);

%x(n+1) = x(n) + dt/6 * (RK-W1 + 2*RK-W2 + 2*RK-W3 + RK-W4);
end


%////////////////////////// differential equation //////////////////////////
%--------- (8.2) ---------
function [ddeltadt] = RK_delta(w,dt,f,n)
	w0 = 2*pi*f;
	ddeltadt = w(n,:)./w0 - 1;
end
%--------- (8.2) ---------

%--------- (8.3) ---------
function [dwdt] = RK_w(w,dt,n,f,H,Pm,Pe,D,numG)
  dwdt = zeros(1,numG);
  w0 = 2*pi*f;
  dwdt = (w0./2./H.*(w0./w(n,:).*Pm(n,:) - w0./w(n,:).*Pe(n,:) - ...
      D.*(w(n,:)./w0 - 1))) ;
end
%--------- (8.3) ---------

%--------- (8.4) ---------
function [deqdt] = RK_eq(numG,Td,ef,xd,xdd,xddd,xl,Kd,eqq,eq,w,id,n)
  deqdt = zeros(1,numG);
  deqdt = 1./Td.*(ef(n,:) + (xd - xdd).*(xdd - xddd)./((xdd - ...
      xl).^2).*Kd.*eqq(n,:) + ...
      (1+(xd-xdd).*(xdd-xddd)./((xdd-xl).^2)).*eq(n,:) + ...
      w(n,:).*(xd-xdd).*(xddd-xl)./(xdd-xl).*id(n,:));
end
%--------- (8.4) ---------

%--------- (8.5) ---------
function [deqqdt] = RK_eqq(numG,Tdd,Kd,eq,eqq,w,xdd,xl,id,n)
  deqqdt = zeros(1,numG);
  deqqdt = ...
      -1./Tdd./Kd.*(Kd.*eqq(n,:)-eq(n,:)+w(n,:).*(xdd-xl).*id(n,:));
end
%--------- (8.5) ---------

%--------- (8.6) ---------
function [deddt] = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd,w,iq,n,ed)
  deddt = zeros(1,numG);
  deddt = -1./Tq.*(-(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2).*Kq.*edd(n,:)...
      +(1+(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2)).*ed(n,:) ...
      -w(n,:).*(xq-xqq).*(xqqq-xl)./(xqq-xl).*iq(n,:));
end
%--------- (8.6) ---------
%--------- (8.7) ---------
function [dedddt] = RK_edd(numG,Tqq,Kq,edd,ed,w,xqq,xl,iq,n)
  dedddt = zeros(1,numG);
  dedddt = -1./Tqq./Kq.*(Kq.*edd(n,:)-ed(n,:)+w(n,:).*(xqq-xl).*iq(n,:));
end
%--------- (8.7) ---------

%--------- (8.115) ---------
function [defdt] = RK_ef(TA,KA,vd,vq,ef0,V,n,ef)
  defdt = -1./TA.*(ef(n,:)+KA.*sqrt(vd(n,:).^2+vq(n,:).^2)...
      -(ef0+KA.*V(1:3))); %V hanyousei
end
%--------- (8.115) ---------

%--------- (8.116) ---------
function [dPmdt] = RK_Pm(Pm,w,numG,TG,KG,w0,n)
  dPmdt = zeros(1,numG);
  dPmdt = -1./TG .* (Pm(n,:) + KG ./ w0 .* w(n,:) - (Pm(1,:) + KG));
end
%--------- (8.116) ---------
%////////////////////////// differential equation //////////////////////////


function [id] = RK_id(IG,delta)
  id = real(IG.*exp(j*(pi/2-delta)));
end

function [IG] = RK_IG()
	IGD = Yg*dEGDQ;
	IG = zeros(1,numG);
	for k = 1:numG
		[tmpTHEATA IG(k)] = cart2pol(IGD(2*k-1),IGD(2*k));
	end
end

function [EGDQ] = RK_EGDQ()
	EGD = zeros(numG,1);
	EGQ = zeros(numG,1);
	EGD = egd .* sin(delta) + egq .* cos(delta);
	EGQ = -egd .* cos(delta) + egq .* sin(delta);
	for k = 1:numG
		EGDQ(2*k-1) = EGD(k);
		EGDQ(2*k) = EGQ(k);
	end
end

function [egd,egq] = RK_egdq(Kq,Kd,xqq,xqqq,xdd,xddd,xl,edd,eqq,ed,eq)
	egd = Kq.*(xqq-xqqq)./(xqq-xl).*edd + (xqqq-xl)./(xqq-xl).*ed;
	egq = Kd.*(xdd-xddd)./(xdd-xl).*eqq + (xddd-xl)./(xdd-xl).*eq;
end

function [Pe] = RK_Pe(egd,id,egq,iq,w0,xddd,xqqq,w)
	Pe = (egd.*id+egq.*iq-w0.*(xddd-xqqq).*id.*iq).*w./w0;
end


function [] = Runge_Kutta(P,numG,Pe,H,D,TG,KG,Td,xd,xdd,xddd,xl,id,Kd,Kd,Kq);


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


for k = 1:numG
  w(1,k) = w0;
  Pm(1,k) = P(k);
end

n = 1;

%--------- k1 ---------
w_k1 = RK_w(w(n,:),dt,n,f,H,Pm,Pe,D,numG)
Pm_k1 = RK_Pm(Pm(n,:),w(n,:),numG,TG,KG,w0,n) 
eq_k1 = RK_eq(numG,Td,ef(n,:),xd,xdd,xddd,xl,...
    Kd,eqq,eq,w(n,:),id(n,:),n)
eqq_k1 = RK_eqq(numG,Tdd,Kd,eq(n,:),w(n,:),xdd,xl,id(n,:)
%--------- k1 ---------
%--------- k2 ---------
wk_2 = RK_w(w(n,:) + dt .* w_k1./2 ,dt,n,f,H,Pm(n,:) + dt .* ...
    Pm_k1./2,Pe,D,numG)
Pm_k2 = RK_Pm(Pm(n,:)+dt.*Pm_k1./2,w(n,:)+dt.*w_k1./2,numG,TG,KG,w0,n)
eq_k2 = RK_eq(numG,Td,ef(n,:)+ef_k1.*dt./2,xd,xdd,xddd,xl,Kd,...
    eqq(n,:)+eqq_k1.*dt./2,eq(n,:)+eq_k1.*dt./2,w(n,:)+w_k1.*dt./2,id(n,:)+id_k1.*dt./2,n)
eqq_k2 = RK_eqq(numG,Tdd,Kd,eq(n,:)+eqq_k1.*dt./2,w(n,:)+w_k1.*dt./2,xdd,xl,id(n,:)+
%--------- k2 ---------
%k3 = RK_w(x(n) + dt/2 * k2);
%k4 = RK_w(x(n) + dt * k3);

%x(n+1) = x(n) + dt/6 * (RK-W1 + 2*RK-W2 + 2*RK-W3 + RK-W4);
end

function [dw] = RK_w(w,dt,n,f,H,Pm,Pe,D,numG)
  w
  dw = zeros(1,numG);
  w0 = 2*pi*f;
  dw = (w0./2./H.*(w(n,:)./w0.*Pm(n,:) - w(n,:)./w0.*Pe(n,:) - ...
      D.*(w(n,:)./w0 - 1))) ;
end

function [dPm] = RK_Pm(Pm,w,numG,TG,KG,w0,n)
  dPm = zeros(1,numG);
  dPm = -1./TG .* (Pm(n,:) + KG ./ w0 .* w(n,:) - (Pm(1,:) + KG));
end

function [deq] = RK_eq(numG,Td,ef,xd,xdd,xddd,xl,Kd,eqq,eq,w,id,n)
  deq = zeros(1,numG);
  deq = 1./Td.*(ef(n,:) + (xd - xdd).*(xdd - xddd)./((xdd - ...
      xl).^2).*Kd.*eqq(n,:) + ...
      (1+(xd-xdd).*(xdd-xddd)./((xdd-xl).^2)).*eq(n,:) + ...
      w(n,:).*(xd-xdd).*(xddd-xl)./(xdd-xl).*id(n,:));
end

function [deqq] = RK_eqq(numG,Tdd,Kd,eq,w,xdd,xl,id)
  deqq = zeros(1,numG);
  deqq = ...
      -1./Tdd./Kd.*(Kd.*eqq(n,:)-eq(n,:)+w(n,:).*(xdd-xl).*id(n,:));
end

function [ded] = RK_ed()
  ded = zeros(1,numG);
  ded = -1./Tq.*(-(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2).*Kq.*edd(n,:)...
      +(1+(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2)).*ed(n,:) ...
      -w(n,:).*(xq-xqq).*(xqqq-xl)./(xqq-xl).*iq(n,:));
end

function [dedd] = RK_edd()
  dedd = zeros(1,numG);
  dedd = -1./Tqq./Kq.*(Kq.*edd(n,:)-ed(n,:)+w(n,:).*(xqq-xl).*iq(n,:));
end

function [def] = RK_ef()
  def = -1./TA.*(ef(n,:)+KA.*sqrt(vd(n,:).^2+vq(n,:).^2)...
      -(ef0+KA.*V(1:3))); %V hanyousei
end

function [did] = RK_id(IG,delta)
  did = real(IG.*exp(j*(pi/2-delta)));
end

function [dIG] = RK_IG()
	IGD = Yg*dEGDQ;
	dIG = zeros(1,numG);
	for k = 1:numG
		[tmpTHEATA dIG(k)] = cart2pol(IGD(2*k-1),IGD(2*k));
	end
end

function [dEGDQ] = RK_EGDQ()
	dEGD = zeros(numG,1);
	dEGQ = zeros(numG,1);
	dEGD = egd .* sin(delta) + egq .* cos(delta);
	dEGQ = -egd .* cos(delta) + egq .* sin(delta);
	for k = 1:numG
		dEGDQ(2*k-1) = dEGD(k);
		dEGDQ(2*k) = dEGQ(k);
	end
end

function [degd,degq] = RK_degdq(Kq,Kd,xqq,xqqq,xdd,xddd,xl,edd,eqq,ed,eq)
	degd = Kq.*(xqq-xqqq)./(xqq-xl).*edd + (xqqq-xl)./(xqq-xl).*ed;
	degq = Kd.*(xdd-xddd)./(xdd-xl).*eqq + (xddd-xl)./(xdd-xl).*eq;
end

function [dPe] = RK_Pe(egd,id,egq,iq,w0,xddd,xqqq,w)
	dPe = (egd.*id+egq.*iq-w0.*(xddd-xqqq).*id.*iq).*w./w0;
end


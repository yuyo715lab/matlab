function [] = Runge_Kutta(P,numG,Pe);


%t(n+1) = t(n) + dt;
max = 200;
f = 50;

dt = 0.01;
w0 = 2*pi*f;
w = zeros(max,numG);
Pm = zeros(max,numG);

for k = 1:numG
  w(1,k) = w0;
  Pm(1,k) = P(k);
end

n = 1;

wk1 = RK_w(w(n,:),dt,n,f,H,Pm,Pe,D,numG)
Pmk1 = RK_Pm(Pm(n,:),w(n,:),numG,TG,KG,w0,n) 

wk2 = RK_w(w(n,:) + dt .* wk1./2 ,dt,n,f,H,Pm(n,:) + dt .* ...
    Pmk1./2,Pe,D,numG)
Pmk2 = RK_Pm(Pm(n,:)+dt.*Pmk1./2,w(n,:)+dt.*wk1./2,numG,TG,KG,w0,n)
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

function [deq] = RK_eq()
  deq = zeros(1,numG);
  deq = 1./Td.*(def(n,:) + (xd - xdd).*(xdd - xddd)./((xdd - ...
      xl).^2).*Kd.*eqq(n,:) + ...
      (1+(xd-xdd).*(xdd-xddd)./((xdd-xl).^2)).*eq(n,:) + ...
      w(n,:).*(xd-xdd).*(xddd-xl)./(xdd-xl).*id(n,:));
end

function [deqq] = RK_eqq()
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

  
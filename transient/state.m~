function [Idq,Vdq,Edq,deltaEq,id,iq,vd,vq,ef0,Pe,eq,eqq,ed,edd,egd,egq,Egd,Egq] = ...
  state(N,numG,GorL,RHO,THEATA,Yprime,Glabel,xd,xq,Rg,Kd,Kq,xdd,xl,xqq,xqqq,xddd)

  
%--------- Id Iq ---------
Vdq = zeros(numG*2,1);
for k = 1:(N - nnz(GorL))
    Vdq(2*k-1) = RHO(Glabel(k)) * cos(THEATA(Glabel(k)));
    Vdq(2*k) = RHO(Glabel(k)) * sin(THEATA(Glabel(k)));
end
Idq = zeros(numG*2,1);
Idq = Yprime * Vdq;
%--------- Id Iq ---------

%--------- Ed Eq ---------
Edq = zeros(numG*2,1);
for k = 1:(N - nnz(GorL))
    Edq(2*k-1) = Vdq(2*k-1) + Rg(k)*Idq(2*k-1) - xq(k) * Idq(2*k);
    Edq(2*k) = Vdq(2*k) + Rg(k)*Idq(2*k) + xq(k) * Idq(2*k-1);
end

%--------- deltaEq ---------
deltaEq = zeros(1,numG);
for k = 1:numG
  deltaEq(k) = cart2pol(Edq(2*k-1),Edq(2*k));
end
%--------- deltaEq ---------

%--------- id iq ---------
id = zeros(1,numG);
iq = zeros(1,numG);
vd = zeros(1,numG);
vq = zeros(1,numG);
for k = 1:numG
  id(k) = real((Idq(2*k-1)+Idq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
  iq(k) = imag((Idq(2*k-1)+Idq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
  vd(k) = real((Vdq(2*k-1)+Vdq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
  vq(k) = imag((Vdq(2*k-1)+Vdq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
end
%--------- id iq ---------

%--------- ef0 ---------
ef0 = zeros(1,numG);
for k = 1:numG
  ef0(k) = vq(k) + xd(k)*id(k) + Rg(k)*iq(k);
end
%--------- ef0 ---------

%--------- Pe ---------
Pe = zeros(1,numG);
Pe = vd .* id + vq .* iq + Rg .* (id.^2 + iq.^2);
%--------- Pe ---------


%--------- eq eqq ed edd --------
eq = zeros(1,numG);
eqq = zeros(1,numG);
ed = zeros(1,numG);
edd = zeros(1,numG);
eq = vq + xdd .* id + Rg .* iq;
eqq = 1 ./ Kd .* (eq - (xdd - xl) .* id);
ed = (xq - xqq) .* iq;
edd = 1 ./ Kq .* (xq - xl) .* iq;
%--------- eq eqq ed edd ---------

%--------- egd egq Egd Egq ---------
egd = zeros(1,numG);
egq = zeros(1,numG);
Egd = zeros(1,numG);
Egq = zeros(1,numG);
egd = Kq .* (xqq - xqqq) ./ (xqq - xl) .* edd + (xqqq - xl) ./ (xqq - ...
    xl) .* ed;
egq = Kd .* (xdd - xddd) ./ (xdd - xl) .* eqq + (xddd - xl) ./ (xdd - ...
    xl) .* eq;
Egd = egd .* sin(deltaEq) + egq .* cos(deltaEq);
Egq = -egd .* cos(deltaEq) + egq .* sin(deltaEq);
%--------- egd egd Egd Egd ---------


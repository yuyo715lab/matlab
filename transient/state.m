function [Idq,Vdq,Edq,deltaEq,id,iq,vd,vq,ef0] = state(N,numG,GorL,RHO,THEATA,Yprime,Glabel,xd,xq,Rg)

  
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

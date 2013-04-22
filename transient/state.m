function [Idq,Edq] = state(N,numG,GorL,RHO,THEATA,Yprime,Glabel,R,xq)

  
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
    Edq(2*k-1) = Vdq(2*k-1) + R(k)*Idq(2*k-1) - xq(k) * Idq(2*k);
    Edq(2*k) = Vdq(2*k) + R(k)*Idq(2*k) + xq(k) * Idq(2*k-1);
end

%--------- deltaEq ---------
deltaEq = zeros(1,numG);
for k = 1:numG
  deltaEq(k) = cart2pol(Edq(2*k-1),Edq(2*k));
end
%--------- deltaEq ---------
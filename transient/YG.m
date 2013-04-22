function [Yg] = YG(numG,xddd,xqqq,Rg,deltaEq,YprimeEF)


Zprime = zeros(2*numG);
for k = 1:numG
  Zprime(2*k-1:2*k,2*k-1:2*k) = [Rg(k) + (xddd(k) - ...
	xqqq(k))*sin(deltaEq(k))*cos(deltaEq(k)) ...
	-(xddd(k)*cos(deltaEq(k))^2 + xqqq(k)*sin(deltaEq(k))^2)
    xddd(k)*sin(deltaEq(k))^2 + xqqq(k)*cos(deltaEq(k))^2 Rg(k) - ...
	(xddd(k) - xqqq(k))*sin(deltaEq(k))*cos(deltaEq(k)) ];
end

Yg = zeros(2*numG);
Yg = inv(eye(2*numG) + YprimeEF*Zprime)*YprimeEF
    
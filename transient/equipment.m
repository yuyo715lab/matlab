function [xd,xdd,xddd,xq,xqq,xqqq,xl,Td,Tdd,Tq,Tqq,Rg,KG,TG,KA,TA,D,H,Kd,Kq] ...
      = equipment(numG)
  
  xd = zeros(1,numG);
  xdd = zeros(1,numG);
  xddd = zeros(1,numG);
  xq = zeros(1,numG);
  xqq = zeros(1,numG);
  xqqq = zeros(1,numG);
  xl = zeros(1,numG);
  Td = zeros(1,numG);
  Tdd = zeros(1,numG);
  Tq = zeros(1,numG);
  Tqq = zeros(1,numG);
  Rg = zeros(1,numG);
  KG = zeros(1,numG);
  TG = zeros(1,numG);
  KA = zeros(1,numG);
  TA = zeros(1,numG);
  D = zeros(1,numG);
  H = zeros(1,numG);
  
  xd = [1.569 1.651 1.220];
  xdd = [0.324 0.232 0.174];
  xddd = [0.249 0.171 0.134];
  xq = [1.548 1.590 1.160];
  xqq = [0.918 0.380 0.250];
  xqqq = [0.248 0.171 0.134];
  xl = [0.204 0.102 0.078];
  Td = [5.14 5.9 8.97];
  Tdd = [0.0437 0.033 0.033];
  Tq = [0.5 0.535 1.5];
  Tqq = [0.07 0.078 0.141];
  H = [50 9 6];
  D = [2.0 2.0 2.0];
  TG = [2.0 2.0 2.0];
  KG = [20.0 20.0 20.0];
  KA = [25.0 25.0 25.0];
  TA = [0.5 0.5 0.5];
  
  Kd = zeros(1,numG);
  Kq = zeros(1,numG);
  Kd = ones(1,numG) + ((xdd - xl) .* (xddd - xl)) ./ ((xdd - xddd) .* ...
      (xd - xl));
  Kq = ones(1,numG) + ((xqq - xl) .* (xqqq - xl)) ./ ((xqq - xqqq) .* ...
      (xq - xl));
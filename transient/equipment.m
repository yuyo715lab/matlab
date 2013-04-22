function [xd,xdd,xddd,xq,xqq,xqqq,xl,Td,Tdd,Tq,Tqq,R,Kg,Tg,Ka,Ta,D,H] ...
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
  R = zeros(1,numG);
  Kg = zeros(1,numG);
  Tg = zeros(1,numG);
  Ka = zeros(1,numG);
  Ta = zeros(1,numG);
  D = zeros(1,numG);
  H = zeros(1,numG);
  
  xd = [1.569 1.651 1.220];
  xdd = [0.324 0.232 0.174];
  xddd = [0.249 0.171 0.134];
  xq = [1.548 1.590 1.160];
  xqq = [0.918 0.380 0.250];
  xqqq = [0.248 0.171 0.134];
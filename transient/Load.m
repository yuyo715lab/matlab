%function Y_ = Load(N,Y,Ps,Qs,PQorPV)

clear all

%--------- power flow calculation ---------
[N,Ps,Qs,PQorPV,P,Q,RHO,THEATA,Y,GorL] = main();
%--------- power flow calculation ---------


%--------- constant impedance load ---------
load_check = zeros(1,N);
for k = 1:N
  if (Ps(k) ~= 0) && (Qs(k) ~= 0 )&& (PQorPV(k) == 1)
    load_check(k) = 1;
  end
end

load_node = zeros(numel(find(load_check)));
load_node = find(load_check);
load_imp = zeros(1,N);

for k = 1:N
  if load_check(k) ~= 0
    load_imp(k) = (-P(k) + Q(k)*i)/RHO(k)^2;
  end
end

Y_ = Y;

for k = 1:numel(load_node)
  Y_(load_node(k),load_node(k)) = Y_(load_node(k),load_node(k)) + load_imp(load_node(k));
end
%--------- constant impedance load ---------

Lnode = find(GorL);
Gnode = find(GorL == 0);
numG = numel(Gnode);
numL = numel(Lnode);
Ynn = zeros(numG);
Ynr = zeros(numG,numL);
Yrr = zeros(numL);
for k = 1:numG
  for m = 1:numG
    Ynn(k,m) = Y_(Gnode(k),Gnode(m));
  end
  for m = 1:numL
    Ynr(k,m) = Y_(Gnode(k),Lnode(m));
  end
end

for k = 1:numL
  for m = 1:numL
    Yrr(k,m) = Y_(Lnode(k),Lnode(m));
  end
  for m = 1:numG
    Yrn(k,m) = Y_(Lnode(k),Gnode(m));
  end
end

Ycont = zeros(numG);
Ycont = Ynn - Ynr*inv(Yrr)*Yrn

Yprime = zeros(2*numG);

for k = 1:numG
  for m = 1:numG
    Yprime(2*k-1:2*k,2*m-1:2*m) = [real(Ycont(k,m)) ...
	  -imag(Ycont(k,m));imag(Ycont(k,m)) real(Ycont(k,m))];
  end
end

Yprime

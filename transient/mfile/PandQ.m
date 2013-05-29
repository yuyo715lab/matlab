function [dP,dQ,dV,P,Q,V] = PandQ(N,Y,PQ,Ps,Qs,Vs,e,f)

Y_ = Y'; %アドミタンス行列の共役転置

for k = 1:N
  for m = 1:N
    PQ(k) = PQ(k) + Y_(k,m) * (e(m) - f(m)*i) * (e(k) + f(k)*i);
  end
end

for k = 1:N
  P(k) = real(PQ(k));
  Q(k) = imag(PQ(k));
end

for k = 1:N
  dP(k) = Ps(k) - P(k);
  dQ(k) = Qs(k) - Q(k);
end

for k = 1:N
  V(k) = e(k)^2 + f(k)^2;
  dV(k) = Vs(k)^2 - V(k);
end

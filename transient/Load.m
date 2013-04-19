%function Y_ = Load(N,Y,Ps,Qs,PQorPV)

load_check = zeros(1,N);
for k = 1:N
  if (Ps(k) ~= 0) && (Qs(k) ~= 0 )&& (PQorPV(k) == 1)
    load_check(k) = 1;
  end
end

load_node = zeros(numel(find(load_check)));
load_node = find(load_check);
load_imp = zeros(N);

for k = 1:N
  if load_check(k) ^= 0
    load_imp(k) = (Ps(k) - Qs(k)*i)/V


function defful = delta(N,PQorPV,NonRef,dP,dQ,dV,Jacobi,e,f)

dPQ = zeros(1,2*(N-1));

for k = 1:N-1
  if PQorPV(NonRef(k)) == 1
    dPQ(2*k-1:2*k) = [dP(NonRef(k)) dQ(NonRef(k))];
  else
    dPQ(2*k-1:2*k) = [dP(NonRef(k)) dV(NonRef(k))];
  end
end


def = zeros(2*(N-1),1);

def = inv(Jacobi) * dPQ.';

defful = zeros(2*N);
m = 1;
for k = 1:N
  if PQorPV(k) == 0
    defful(2*k-1:2*k,1) = 0;
  else
    defful(2*k-1:2*k,1) = def(m:m+1,1);
    m = m + 2;
  end
end






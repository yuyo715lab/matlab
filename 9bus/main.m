%/////// 9-bus power flow calculation
%/////// TOMOHIRO ADACHI
%/////// 2013/04


clear all


%============= Define.m ================
[N,Ref,PQorPV,NonRef,R,Tr,e,f,Vs,V,dV,Ps,Qs,PQ] = Define();
%============= Define.m ================

%============= admittance.m ================
Y = admittance(N,R,Tr);
%============= admittance.m ================

%============= PandQ.m ================
[dP,dQ,dV,P,Q,V] = PandQ(N,Y,PQ,Ps,Qs,Vs,e,f);
%============= PandQ.m ================

%============= jacobi.m ================
Jacobi = jacobi(N,Y,e,f,PQorPV,NonRef);
%============= jacobi.m ================

%============= delta.m ================
def = delta(N,PQorPV,NonRef,dP,dQ,dV,Jacobi,e,f);
%============= delta.m ================


n = 1 % the number of cycle

%------------ polar coordinate system ------------
THEATA = zeros(1,N-1); %theta
RHO = zeros(1,N-1); %r
for k = 1:N-1
  [THEATA(k),RHO(k)] = cart2pol(e(NonRef(k)) + def(2*k-1),(f(NonRef(k)) + def(2*k)));
end
%------------ polar coordinate system ------------


%/////// Display /////////
THEATA
RHO
P
Q
%/////// Display /////////

b_rho = zeros(1,N-1);
a_rho = RHO;
d_rho = zeros(1,N-1);
for k = 1:N-1
  d_rho(k) = abs(a_rho(k) - b_rho(k)); %|V - V|
end

while max(d_rho) > 10^(-7) %convergence conditiongit
  n = n + 1
  for k = 1:N-1
    e(NonRef(k)) = e(NonRef(k)) + def(2*k-1);
    f(NonRef(k)) = f(NonRef(k)) + def(2*k);
  end
  
  %============= PandQ.m ================
  [dP,dQ,dV,P,Q,V] = PandQ(N,Y,PQ,Ps,Qs,Vs,e,f);
  %============= PandQ.m ================
  
  %============= jacobi.m ================
  Jacobi = jacobi(N,Y,e,f,PQorPV,NonRef);
  %============= jacobi.m ================
  
  %============= delta.m ================
  def = delta(N,PQorPV,NonRef,dP,dQ,dV,Jacobi,e,f);
  %============= delta.m ================
  
  %----------- polar coordinate system ---------
  for k = 1:N-1
    [THEATA(k),RHO(k)] = cart2pol(e(NonRef(k)) + def(2*k-1),(f(NonRef(k)) + def(2*k)));
  end
  %----------- polar coordinate system ---------
  
  
  b_rho = a_rho;
  a_rho = RHO;
  for k = 1:N-1
    d_rho(k) = abs(a_rho(k) - b_rho(k));
  end

%/////// Display /////////
THEATA
RHO
P
Q
%/////// Display /////////



end


    

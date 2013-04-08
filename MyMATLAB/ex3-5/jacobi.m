function [Jacobi] = jacobi(N,Y,e,f,PQorPV,NonRef)

G = real(Y);
B = imag(Y);
a = zeros(1,N);
b = zeros(1,N);
Pe = zeros(N,N);
Qf = zeros(N,N);
Pf = zeros(N,N);
Qe = zeros(N,N);
Ve = zeros(N,N);
Vf = zeros(N,N);
Jacobi = zeros(2*(N-1),2*(N-1));  


for k = 1:N
  for m = 1:N
    a(k) = a(k) + G(k,m) * e(m) - B(k,m) * f(m);
    b(k) = b(k) + G(k,m) * f(m) + B(k,m) * e(m);
  end
end


for k = 1:N
  for m = 1:N
    if k ~= m
      Pe(k,m) = G(k,m) * e(k) + B(k,m) * f(k);
      Qf(k,m) = -Pe(k,m);
      Pf(k,m) = -B(k,m) * e(k) + G(k,m) * f(k);
      Qe(k,m) = Pf(k,m);
      Ve(k,m) = 0;
      Vf(k,m) = 0;
    else
      Pe(k,m) = a(k) + G(k,m) * e(k) + B(k,m) * f(k);
      Qf(k,m) = a(k) - G(k,m) * e(k) - B(k,m) * f(k);
      Pf(k,m) = b(k) - B(k,m) * e(k) + G(k,m) * f(k);
      Qe(k,m) = -b(k) - B(k,m) * e(k) + G(k,m) * f(k);
      Ve(k,m) = 2 * e(k);
      Vf(k,m) = 2 * f(k);
    end
  end
end



Jrow = 1;
Jline = 1;
%%NonRef = find(PQorPV);

for k = 1:N
  if PQorPV(k) == 1
    for m = 1:length(NonRef)
      Jacobi(Jrow:Jrow+1,Jline:Jline+1) ... 
	  = [Pe(k,NonRef(m)) Pf(k,NonRef(m))
	Qe(k,NonRef(m)) Qf(k,NonRef(m))];
      Jline = Jline + 2;
    end
    Jline = 1;
    Jrow = Jrow + 2;
  else 
    if PQorPV(k) == 2
      for m = 1:length(NonRef)
	Jacobi(Jrow:Jrow+1, Jline:Jline+1) = [Pe(k,NonRef(m)) ...
	      Pf(k,NonRef(m));Ve(k,NonRef(m)) Vf(k,NonRef(m))];
	Jline = Jline + 2;
      end
      Jline = 1;
      Jrow = Jrow + 2;
    end
  end
end


Jrow = 1;
Jline = 1;

    
function Out = admittance(N,R,Tr)

Y = zeros(N,N);

%%%%%%%%% Direct admittance %%%%%%%%%
for k = 1:N
  for m = 1:N
    if R(k,m) ~= 0
      Y(k,k) = Y(k,k) + 1/R(k,m);
    end
  end
end
%%%%%%%%% Direct admittance %%%%%%%%%


%%%%%%%%% Transfer admittance %%%%%%%%%
for k = 1:N
  for m = k:N
    if k ~= m
      if R(k,m) ~= 0
	Y(k,m) = -1/R(k,m);
      end
    end
    Y(m,k) = Y(k,m);
  end
end
%%%%%%%%% Transfer admittance %%%%%%%%%

%%%%%%%%% Off norminal turn rario %%%%%%%%%
for k = 1:N
  for m = k:N
    if Tr(k,m) ~= 0
      if R(k,m) ~= 0 && R(m,m) ~= 0
	Y(k,k) = Y(k,k) + (Tr(k,m)^2 - 1) * (1/R(k,m) + 1/R(m,m));
	Y(k,m) = Y(k,m) - (Tr(k,m) - 1) * (1/R(k,m) + 1/R(m,m));
      else 
	if R(k,m) == 0 && R(m,m) ~= 0
	  Y(k,k) = Y(k,k) + (Tr(k,m)^2 - 1) * 1/R(m,m);
	  Y(k,m) = Y(k,m) - (Tr(k,m) - 1) * 1/R(m,m);
	else
	  if R(k,m) ~= 0 && R(m,m) == 0
	    Y(k,k) = Y(k,k) + (Tr(k,m)^2 - 1) * 1/R(k,m);
	    Y(k,m) = Y(k,m) - (Tr(k,m) - 1) * 1/R(k,m);
	  end
	end
      end 
    end
    Y(m,k) = Y(k,m);
  end
end
%%%%%%%%% Off norminal turn rario %%%%%%%%%


Out = Y;

      
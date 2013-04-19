n = numel(XY)/4 - 1;

dvdp = zeros(n,2);

for k = 1:n
  dvdp(k,1) = XY(k,1);
  dvdp(k,2) = (XY(k+1,2)-XY(k,2))/(XY(k+1,1)-XY(k,1));
end

scatter(dvdp(:,1),dvdp(:,2))
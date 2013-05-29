r = 1;
for k = 2:max
	if (delta_for_plot(k,2) - delta_for_plot(k-1,2)) > 0 && ...
			 (delta_for_plot(k,2) - delta_for_plot(k+1,2)) > 0
		maxTime(r,1) = k*dt;
		maxTime(r,2) = delta_for_plot(k,2)-delta_for_plot(1,2);
		r = r+1;
	end
end
r = 1;
bar(maxTime(:,1),maxTime(:,2),'BarWidth',0.3);

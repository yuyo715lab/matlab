max = round(13/0.00001)+1;
delta_for_plot = zeros(max,6);
delta_for_plot(:,1) = 0:0.00001:13;
delta_for_plot(:,2) = dlmread('./csv/delta_4142_00001.csv',' ',[0 1 max-1 1]);
delta_for_plot(:,3) = dlmread('./csv/delta_4142_00001.csv',' ',[0 2 max-1 2]);
delta_for_plot(:,4) = dlmread('./csv/delta_43_00001.csv',' ',[0 1 max-1 1]);
delta_for_plot(:,5) = dlmread('./csv/delta_44_00001.csv',' ',[0 1 max-1 1]);
delta_for_plot(:,6) = dlmread('./csv/delta_45_00001.csv',' ',[0 1 max-1 1]);
plot(delta_for_plot(:,1),delta_for_plot(:,2),delta_for_plot(:,1),delta_for_plot(:,3),...
	 delta_for_plot(:,1),delta_for_plot(:,4),delta_for_plot(:,1),...
	 delta_for_plot(:,5),delta_for_plot(:,1),delta_for_plot(:,6),'LineWidth',2)
grid on
xlabel('time[sec]')
ylabel('phase difference angle[degree]')
legend('0.41','0.42','0.43','0.44','0.45')

%************** transient calculation
%************** Tomohiro Adachi
%************** 2013/05

clear all
format long;
startT = clock();
startCpuT = cputime;

%././././././--------------- ././././././
%././././././ YOU CAN EDIT  ././././././
dt = 0.01; % sampling time
endTime = 30;
removeAccidentYesNo = 1;
removeTime = 0.2;
EarthFault_from = 5;
EarthFault_to = 7;
OpenYesNo = 0;
OpenTime = 0.5;
plotColor = 3; % ['r' 'g' 'c' 'y' 'm' 'b' 'k'];
Hold = 0; % hold  on > 1   off > 0
step_mail_yesno = 0; % yes > 1 no > 0
finish_mail_yesno = 0; % yes > 1 no > 0
delta_w_v = 0; % 0 delta > 1 frequency > 2 Voltage > 3
csvname = '../csv/test.csv';
%././././././ YOU CAN EDIT ././././././
%././././././---------------././././././


%:::::::::::::::: Don't Edit below!!! ::::::::::::::::::::::::::::: 
%|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
%|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
max = round(endTime/dt);

%--------- power flow calculation ---------
[N,Ps,Qs,PQorPV,P,Q,RHO,THEATA,Y,GorL] = PowerCalc();
%--------- power flow calculation ---------

%--------- Yprime ---------
[YprimeStart,numG,numL] = Yprime(N,Y,Ps,Qs,PQorPV,P,Q,RHO,GorL);
%--------- Yprime ---------

%--------- YprimeEF ---------
[YprimeEF] = YprimeEF(N,Y,Ps,Qs,PQorPV,P,Q,RHO,GorL,EarthFault_from);
%--------- YprimeEF ---------

%--------- YprimeOpen ---------
YprimeOpen = open(EarthFault_from,EarthFault_to,P,Q,RHO,GorL);
%--------- YprimeOpen ---------

%--------- eqipment constant ---------
[xd,xdd,xddd,xq,xqq,xqqq,xl,Td,Tdd,Tq,Tqq,Rg,KG,TG,KA,TA,D,H,Kd,Kq] ...
	 = equipment(numG);
%--------- eqipment constant ---------

%--------- state ---------
Glabel = zeros(1,numG);
Glabel = find(GorL == 0);
[Idq,Vdq,Edq,deltaEq,id,iq,vd0,vq0,ef0,Pe,eq,eqq,ed,edd,egd,egq,Egd,Egq] = ...
	 state(N,numG,GorL,RHO,THEATA,YprimeStart,Glabel,xd,xq,Rg,Kd,Kq,xdd,xl,xqq,xqqq,xddd,max);
%--------- state ---------

%--------- Yg ---------
[Yg] = YG(numG,xddd,xqqq,Rg,deltaEq,YprimeEF);
%--------- Yg ---------

%--------- id iq ---------
Egdq = zeros(numG*2,1);
for k = 1:numG
	Egdq(2*k-1) = Egd(k);
	Egdq(2*k) = Egq(k);
end
Idq = Yg * Egdq;
Vdq = inv(YprimeEF) * Idq;
for k = 1:numG
	id(k) = real((Idq(2*k-1)+Idq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
	iq(k) = imag((Idq(2*k-1)+Idq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
	vd(k) = real((Vdq(2*k-1)+Vdq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
	vq(k) = imag((Vdq(2*k-1)+Vdq(2*k)*i)*exp(i*(pi/2 - deltaEq(k))));
end
%--------- id iq ---------

%--------- Pe ---------
Pe(1,:) = vd .* id + vq .* iq + Rg .* (id.^2 + iq.^2);
%--------- Pe ---------

%!!!!!!!!!!!!!!!!!!!!  Runge Kutta !!!!!!!!!!!!!!!!!!!!!
delta_for_plot = zeros(max+1,2);
w_for_plot = zeros(max+1,2);
v_for_plot = zeros(max+1,2);
now_step = 2;

[delta_for_plot,w_for_plot,v_for_plot] =...
	 Runge_Kutta3(P,numG,Pe,H,D,TG,KG,Td,Tdd,Tq,xd,xdd,xddd,xl,id,Kd,Kq,vd,vq,KA,TA,...
	 xq,xqq,xqqq,iq,Tqq,ef0,deltaEq,eq,eqq,ed,edd,vd0,vq0,Yg,YprimeEF,max,YprimeStart,...
	 Rg,Glabel,removeTime,delta_for_plot,w_for_plot,v_for_plot,dt,now_step,endTime,step_mail_yesno,...
	 YprimeOpen,OpenTime,OpenYesNo,removeAccidentYesNo);
%!!!!!!!!!!!!!!!!!!!!  Runge Kutta !!!!!!!!!!!!!!!!!!!!!


%============== Plot =================
plot_col = ['r' 'g' 'c' 'y' 'm' 'b' 'k'];
plot_col_hasen = ['r:';'g:';'c:';'y:';'m:';'b:';'k:'];
switch delta_w_v
	case 0
		plot(delta_for_plot(:,1),delta_for_plot(:,now_step),plot_col(plotColor),'LineWidth',2)
	case 1
		plot(w_for_plot(:,1),w_for_plot(:,now_step),plot_col(plotColor),'LineWidth',2)
	case 2
		plot(v_for_plot(:,1),v_for_plot(:,now_step),plot_col(plotColor),'LineWidth',2)
end
xlabel('time[sec]')
switch delta_w_v
	case 0 
		ylabel('Phase Difference Angle[degree]')
	case 1
		ylabel('Frequency[Hz]')
	case 2
		ylabel('Voltage[p.u.]')
end
if removeAccidentYesNo == 1
	label = ['remove ' num2str(removeTime) '[sec]'];
else 
	if OpenYesNo == 1
		label = ['open ' num2str(OpenTime) '[sec]'];
	end
end
legend(label)
grid on
if Hold
	hold on
else
	hold off
end
switch delta_w_v
	case 0
		dlmwrite(csvname,delta_for_plot,' ');
	case 1
		dlmwrite(csvname,w_for_plot,' ');
	case 2
		dlmwrite(csvname,v_for_plot,' ');
end
%!sudo chmod a+w ./csv/delta_43_00001.csv
%============== Plot =================

%TTTTTTTTTTTTT Cal time TTTTTTTTTTTTTTT
ntime=cputime-startCpuT;
nhour = floor(ntime/60/60);
nmin = floor((ntime-nhour*3600)/60);
nsec = ntime-nhour*3600-nmin*60;
disp(sprintf('%s%s', 'start time:',datestr(startT,31)));
disp(sprintf('%s%s', 'finish time:',datestr(clock,31)));
disp(sprintf('%s%d%s%02d%s%04.1f%s', ...
	 'calculation time:',nhour,'h',nmin,'m',nsec,'s'));
%TTTTTTTTTTTTT Cal time TTTTTTTTTTTTTTT

if finish_mail_yesno 
	finish_mail()
end

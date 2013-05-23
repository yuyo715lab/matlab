clear all
format long;
startT = clock();
startCpuT = cputime;

%matlabpool 4

dt = 0.01; % sampling time
endTime = 30;
removeAccidentYesNo = 0;
eft_min = 0.2;
eft_max = 0.2;
eft_step = 0.1;
EarthFault_from = 5;
EarthFault_to = 7;
OpenYesNo = 1;
OpenTime = 0.07;
step_mail_yesno = 0; %1 yes 0 no
csvname = './csv/test.csv';

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
%YprimeEF = Yprime; %non EF

Glabel = zeros(1,numG);
Glabel = find(GorL == 0);

%--------- eqipment constant ---------
[xd,xdd,xddd,xq,xqq,xqqq,xl,Td,Tdd,Tq,Tqq,Rg,KG,TG,KA,TA,D,H,Kd,Kq] ...
	 = equipment(numG);
%--------- eqipment constant ---------

%--------- state ---------
[Idq,Vdq,Edq,deltaEq,id,iq,vd0,vq0,ef0,Pe,eq,eqq,ed,edd,egd,egq,Egd,Egq] = ...
	 state(N,numG,GorL,RHO,THEATA,YprimeStart,Glabel,xd,xq,Rg,Kd,Kq,xdd,xl,xqq,xqqq,xddd,max);
%--------- state ---------


%--------- Yg ---------
[Yg] = YG(numG,xddd,xqqq,Rg,deltaEq,YprimeEF);
%--------- Yg ---------

Egdq = zeros(numG*2,1);
for k = 1:numG
	Egdq(2*k-1) = Egd(k);
	Egdq(2*k) = Egq(k);
end

Idq = Yg * Egdq;
Vdq = inv(YprimeEF) * Idq;

%--------- id iq ---------
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


number_of_step = (eft_max > eft_min)*...
	 (round((eft_max - eft_min)/eft_step)+1) ...
	 + (eft_max == eft_min)*1;
delta_for_plot = zeros(max+1,number_of_step+1);
w_for_plot = zeros(max+1,number_of_step+1);
v_for_plot = zeros(max+1,number_of_step+1);
plot_col = ['r' 'g' 'c' 'y' 'm' 'b' 'k'];
plot_col_hasen = ['r:';'g:';'c:';'y:';'m:';'b:';'k:'];
now_step = 2;
label = blanks(5);
for k = eft_min:eft_step:eft_max
	EarthFaultTime = k;
	if mod(EarthFaultTime,0.1) == 0 && mod(EarthFaultTime,1) ~= 0;label = [label;num2str(k),'00'];
	else  
		if mod(EarthFaultTime,1) == 0;label = [label;num2str(k),'.000'];
			else
			label = [label;num2str(k)];end 
			end
			tic
	[delta_for_plot,w_for_plot,v_for_plot] =...
		 Runge_Kutta3(P,numG,Pe,H,D,TG,KG,Td,Tdd,Tq,xd,xdd,xddd,xl,id,Kd,Kq,vd,vq,KA,TA,...
		 xq,xqq,xqqq,iq,Tqq,ef0,deltaEq,eq,eqq,ed,edd,vd0,vq0,Yg,YprimeEF,max,YprimeStart,...
		 Rg,Glabel,EarthFaultTime,delta_for_plot,w_for_plot,v_for_plot,dt,now_step,endTime,step_mail_yesno,...
		 YprimeOpen,OpenTime,OpenYesNo,removeAccidentYesNo);
	
	if now_step < 9
		plot(delta_for_plot(:,1),delta_for_plot(:,now_step),plot_col(now_step-1),'LineWidth',2)
	else
		plot(delta_for_plot(:,1),delta_for_plot(:,now_step),plot_col_hasen(now_step-8,:),'LineWidth',2)
	end
	hold on
	
	now_step = now_step + 1;
end

%PPPPPPPPPPPP Plot PPPPPPPPPPPPPPPPPP
xlabel('time[sec]')
ylabel('phase difference angle[degree]')
trueLabel = label(2:end,:);
legend(trueLabel)
grid on
hold off
%hold on
dlmwrite(csvname,delta_for_plot,' ');

%dlmwrite(csvname,w_for_plot,' ');
%dlmwrite(csvname,v_for_plot,' ');
%!sudo chmod a+w ./csv/delta_43_00001.csv
%PPPPPPPPPPPP Plot PPPPPPPPPPPPPPPPPP


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

%matlabpool close
%finish_mail()

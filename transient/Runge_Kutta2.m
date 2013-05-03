function [] =Runge_Kutta2(P,numG,Pe,H,D,TG,KG,Td,Tdd,Tq,xd,xdd,xddd,xl,id,Kd,Kq,vd,vq,KA,TA,...
		 xq,xqq,xqqq,iq,Tqq,ef0,deltaEq,eq_state,eqq_state,ed_state,edd_state,vd0,vq0,Yg,YprimeEF,...
		 max,Yprime,Rg,Glabel);

	f = 50; % frequency
	dt = 0.01; % sampling time
	w0 = 2*pi*f; % angular velocity

	%$$$$$$$$$ initialization $$$$$$$$$
	w = zeros(max,numG);
	Pm = zeros(max,numG);
	eq = zeros(max,numG);
	ef = zeros(max,numG);
	eqq = zeros(max,numG);
	eqqq = zeros(max,numG);
	ed = zeros(max,numG);
	edd = zeros(max,numG);
	eddd = zeros(max,numG);
	V = zeros(max,numG);
	delta = zeros(max,numG);
	%$$$$$$$$$ initialization $$$$$$$$$

	%((((((((( initial value )))))))))
	delta(1,:) = deltaEq;
	eq(1,:) = eq_state;
	eqq(1,:) = eqq_state;
	ed(1,:) = ed_state;
	edd(1,:) = edd_state;
	ef(1,:) = ef0;
	for k = 1:numG
		w(1,k) = w0;
		Pm0(1,k) = P(Glabel(k));
		Pm(1,k) = P(Glabel(k));
		V(1,k) = sqrt(vd0(k)^2 + vq0(k)^2);
	end
	%((((((((( initial value )))))))))

	%/////////////////////////////////////////////
	%//////////// for loop ///////////////////////
	for n = 1:max-1
		%<<<<<<<<<<<<<< earth fault >>>>>>>>>>>>>>>
		if n == 20 % earth fault occur at n*dt
			[Yg] = YG(numG,xddd,xqqq,Rg,deltaEq,Yprime);
			YprimeEF = Yprime;
		end
		%<<<<<<<<<<<<<< earth fault >>>>>>>>>>>>>>>
		delta_k = zeros(5,numG);
		w_k = zeros(5,numG);
		eq_k = zeros(5,numG);
		eqq_k = zeros(5,numG);
		ed_k = zeros(5,numG);
		edd_k = zeros(5,numG);
		ef_k = zeros(5,numG);
		Pm_k = zeros(5,numG);
		dt_k_tmp = [0;dt/2;dt/2;dt];
		dt_k = dt_k_tmp*ones(1,3);
		id_atk = zeros(5,numG);
		iq_atk = zeros(5,numG);
		vd_atk = zeros(5,numG);
		vq_atk = zeros(5,numG);
		for k = 1:4
			id_atk(k+1,:) = RK_id(delta(n,:)+dt_k(k,:).*delta_k(k,:),Yg,...
				 RK_EGDQ(numG,delta(n,:)+dt_k(k,:).*delta_k(k,:),...
				 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:)),...
				 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt_k(k,:).*eqq_k(k,:),eq(n,:)+dt_k(k,:).*eq_k(k,:))),...
				 numG);
			iq_atk(k+1,:) = RK_iq(delta(n,:)+dt_k(k,:).*delta_k(k,:),Yg,...
				 RK_EGDQ(numG,delta(n,:)+dt_k(k,:).*delta_k(k,:),...
				 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:)),...
				 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt_k(k,:).*eqq_k(k,:),eq(n,:)+dt_k(k,:).*eq_k(k,:))),...
				 numG);
			vd_atk(k+1,:) = RK_vd(delta(n,:)+dt_k(k,:).*delta_k(k,:),Yg,...
				 RK_EGDQ(numG,delta(n,:)+dt_k(k,:).*delta_k(k,:),...
				 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:)),...
				 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt_k(k,:).*eqq_k(k,:),eq(n,:)+dt_k(k,:).*eq_k(k,:))),...
				 numG,YprimeEF);
			vq_atk(k+1,:) = RK_vq(delta(n,:)+dt_k(k,:).*delta_k(k,:),Yg,...
				 RK_EGDQ(numG,delta(n,:)+dt_k(k,:).*delta_k(k,:),...
				 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:)),...
				 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt_k(k,:).*eqq_k(k,:),eq(n,:)+dt_k(k,:).*eq_k(k,:))),...
				 numG,YprimeEF);
			delta_k(k+1,:) = RK_delta(w(n,:)+dt_k(k,:).*delta_k(k,:),dt,f,n);
			w_k(k+1,:) = RK_w(w(n,:)+dt_k(k,:).*w_k(k,:),dt,n,f,H,Pm(n,:)+dt_k(k,:).*Pm_k(k,:),...
				 RK_Pe(RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:)),...
				 id_atk(k+1,:),RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt_k(k,:).*eqq_k(k,:),eq(n,:)+dt_k(k,:).*eq_k(k,:)),...
				 iq_atk(k+1,:),w0,xddd,xqqq,w(n,:)+dt_k(k,:).*w_k(k,:)),D,numG);
			eq_k(k+1,:) = RK_eq(numG,Td,ef(n,:)+dt_k(k,:).*ef_k(k,:),xd,xdd,xddd,xl,Kd,...
				 eqq(n,:)+dt_k(k,:).*eqq_k(k,:),...
				 eq(n,:)+dt_k(k,:).*eq_k(k,:),w(n,:)+dt_k(k,:).*w_k(k,:),id_atk(k+1,:),n,w0);
			eqq_k(k+1,:) =RK_eqq(numG,Tdd,Kd,eq(n,:)+dt_k(k,:).*eq_k(k,:),eqq(n,:)+dt_k(k,:).*eqq_k(k,:),...
				 w(n,:)+dt_k(k,:).*w_k(k,:),xdd,xl,id_atk(k+1,:),n,w0);
			ed_k(k+1,:) = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd(n,:)+dt_k(k,:).*edd_k(k,:),...
				 w(n,:)+dt_k(k,:).*w_k(k,:),...
				 iq_atk(k+1,:),n,ed(n,:)+dt_k(k,:).*ed_k(k,:),w0);
			edd_k(k+1,:)=RK_edd(numG,Tqq,Kq,edd(n,:)+dt_k(k,:).*edd_k(k,:),ed(n,:)+dt_k(k,:).*ed_k(k,:),...
				 w(n,:)+dt_k(k,:).*w_k(k,:),xqq,xl,iq_atk(k+1,:),n,w0);
			ef_k(k+1,:) = RK_ef(TA,KA,vd_atk(k+1,:),vq_atk(k+1,:),ef0,V,n,ef(n,:)+dt_k(k,:).*ef_k(k,:));
			Pm_k(k+1,:)=RK_Pm(Pm(n,:)+dt_k(k,:).*Pm_k(k,:),w(n,:)+dt_k(k,:).*w_k(k,:),numG,TG,KG,w0,n,Pm0);
		end
		if n==1
			delta_k
			w_k
			eq_k
			eqq_k
			ed_k
			edd_k
			ef_k
			Pm_k
		end
		
%{
		%--------- k1 ---------
		delta_k1 = RK_delta(w(n,:),dt,f,n);
		w_k1 = RK_w(w(n,:),dt,n,f,H,Pm(n,:),Pe(n,:),D,numG);
		eq_k1 = RK_eq(numG,Td,ef(n,:),xd,xdd,xddd,xl,...
			 Kd,eqq(n,:),eq(n,:),w(n,:),id(n,:),n,w0);
		eqq_k1 = RK_eqq(numG,Tdd,Kd,eq(n,:),eqq(n,:),w(n,:),xdd,xl,id(n,:),n,w0);
		ed_k1 = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd(n,:),w(n,:),iq(n,:),n,ed(n,:),w0);
		edd_k1 = RK_edd(numG,Tqq,Kq,edd(n,:),ed(n,:),w(n,:),xqq,xl,iq(n,:),n,w0);
		ef_k1 = RK_ef(TA,KA,vd(n,:),vq(n,:),ef0,V,n,ef(n,:));
		Pm_k1 = RK_Pm(Pm(n,:),w(n,:),numG,TG,KG,w0,n,Pm0);
		%--------- k1 ---------
		%--------- k2 ---------
		id_atk2 = RK_id(delta(n,:)+dt./2.*delta_k1,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k1,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k1,ed(n,:)+dt./2.*ed_k1),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k1,eq(n,:)+dt./2.*eq_k1)),...
			 numG);
		iq_atk2 = RK_iq(delta(n,:)+dt./2.*delta_k1,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k1,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k1,ed(n,:)+dt./2.*ed_k1),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k1,eq(n,:)+dt./2.*eq_k1)),...
			 numG);
		vd_atk2 = RK_vd(delta(n,:)+dt./2.*delta_k1,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k1,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k1,ed(n,:)+dt./2.*ed_k1),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k1,eq(n,:)+dt./2.*eq_k1)),...
			 numG,YprimeEF);
		vq_atk2 = RK_vq(delta(n,:)+dt./2.*delta_k1,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k1,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k1,ed(n,:)+dt./2.*ed_k1),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k1,eq(n,:)+dt./2.*eq_k1)),...
			 numG,YprimeEF);
		delta_k2 = RK_delta(w(n,:)+dt./2.*w_k1,dt,f,n);
		w_k2 = RK_w(w(n,:)+dt./2.*w_k1,...
			 dt,n,f,H,...
			 Pm(n,:)+dt./2.*Pm_k1,...
			 RK_Pe(RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k1,ed(n,:)+dt./2.*ed_k1),...
			 id_atk2,...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k1,eq(n,:)+dt./2.*eq_k1),...
			 iq_atk2,w0,xddd,xqqq,w(n,:)+dt./2.*w_k1),D,numG);
		eq_k2 = RK_eq(numG,Td,ef(n,:)+dt./2.*ef_k1,xd,xdd,xddd,xl,Kd,eqq(n,:)+dt./2.*eqq_k1,...
			 eq(n,:)+dt./2.*eq_k1,w(n,:)+dt./2.*w_k1,id_atk2,n,w0);
		eqq_k2 =RK_eqq(numG,Tdd,Kd,eq(n,:)+dt./2.*eq_k1,eqq(n,:)+dt./2.*eqq_k1,w(n,:)+dt./2.*w_k1,xdd,xl,...
			 id_atk2,n,w0);
		ed_k2 = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd(n,:)+dt./2.*edd_k1,w(n,:)+dt./2.*w_k1,iq_atk2,...
			 n,ed(n,:)+dt./2.*ed_k1,w0);
		edd_k2 =RK_edd(numG,Tqq,Kq,edd(n,:)+dt./2.*edd_k1,ed(n,:)+dt./2.*ed_k1,w(n,:)+dt./2.*w_k1,xqq,xl,...
			 iq_atk2,n,w0);
		ef_k2 = RK_ef(TA,KA,vd_atk2,vq_atk2,ef0,V,n,ef(n,:)+dt./2.*ef_k1);
		Pm_k2 = RK_Pm(Pm(n,:)+dt./2.*Pm_k1,w(n,:)+dt./2.*w_k1,numG,TG,KG,w0,n,Pm0);
		%--------- k2 ---------
		%--------- K3 ---------
		id_atk3 = RK_id(delta(n,:)+dt./2.*delta_k2,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k2,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k2,ed(n,:)+dt./2.*ed_k2),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k2,eq(n,:)+dt./2.*eq_k2)),...
			 numG);
		iq_atk3 = RK_iq(delta(n,:)+dt./2.*delta_k2,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k2,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k2,ed(n,:)+dt./2.*ed_k2),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k2,eq(n,:)+dt./2.*eq_k2)),...
			 numG);
		vd_atk3 = RK_vd(delta(n,:)+dt./2.*delta_k2,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k2,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k2,ed(n,:)+dt./2.*ed_k2),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k2,eq(n,:)+dt./2.*eq_k2)),...
			 numG,YprimeEF);
		vq_atk3 = RK_vq(delta(n,:)+dt./2.*delta_k2,Yg,RK_EGDQ(numG,delta(n,:)+dt./2.*delta_k2,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k2,ed(n,:)+dt./2.*ed_k2),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k2,eq(n,:)+dt./2.*eq_k2)),...
			 numG,YprimeEF);
		delta_k3 = RK_delta(w(n,:)+dt./2.*w_k2,dt,f,n);
		w_k3 = RK_w(w(n,:)+dt./2.*w_k2,...
			 dt,n,f,H,...
			 Pm(n,:)+dt./2.*Pm_k2,...
			 RK_Pe(RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt./2.*edd_k2,ed(n,:)+dt./2.*ed_k2),...
			 id_atk3,...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt./2.*eqq_k2,eq(n,:)+dt./2.*eq_k2),...
			 iq_atk3,w0,xddd,xqqq,w(n,:)+dt./2.*w_k2),D,numG);
		eq_k3 = RK_eq(numG,Td,ef(n,:)+dt./2.*ef_k2,xd,xdd,xddd,xl,Kd,eqq(n,:)+dt./2.*eqq_k2,...
			 eq(n,:)+dt./2.*eq_k2,w(n,:)+dt./2.*w_k2,id_atk3,n,w0);
		eqq_k3 =RK_eqq(numG,Tdd,Kd,eq(n,:)+dt./2.*eq_k2,eqq(n,:)+dt./2.*eqq_k2,w(n,:)+dt./2.*w_k2,xdd,xl,...
			 id_atk3,n,w0);
		ed_k3 = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd(n,:)+dt./2.*edd_k2,w(n,:)+dt./2.*w_k2,iq_atk3,...
			 n,ed(n,:)+dt./2.*ed_k2,w0);
		edd_k3 =RK_edd(numG,Tqq,Kq,edd(n,:)+dt./2.*edd_k2,ed(n,:)+dt./2.*ed_k2,w(n,:)+dt./2.*w_k2,xqq,xl,...
			 iq_atk3,n,w0);
		ef_k3 = RK_ef(TA,KA,vd_atk3,vq_atk3,ef0,V,n,ef(n,:)+dt./2.*ef_k2);
		Pm_k3 = RK_Pm(Pm(n,:)+dt./2.*Pm_k2,w(n,:)+dt./2.*w_k2,numG,TG,KG,w0,n,Pm0);
		%--------- k3 ---------
		%--------- K4 ---------
		id_atk4 = RK_id(delta(n,:)+dt.*delta_k3,Yg,RK_EGDQ(numG,delta(n,:)+dt.*delta_k3,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt.*edd_k3,ed(n,:)+dt.*ed_k3),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt.*eqq_k3,eq(n,:)+dt.*eq_k3)),...
			 numG);
		iq_atk4 = RK_iq(delta(n,:)+dt.*delta_k3,Yg,RK_EGDQ(numG,delta(n,:)+dt.*delta_k3,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt.*edd_k3,ed(n,:)+dt.*ed_k3),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt.*eqq_k3,eq(n,:)+dt.*eq_k3)),...
			 numG);
		vd_atk4 = RK_vd(delta(n,:)+dt.*delta_k3,Yg,RK_EGDQ(numG,delta(n,:)+dt.*delta_k3,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt.*edd_k3,ed(n,:)+dt.*ed_k3),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt.*eqq_k3,eq(n,:)+dt.*eq_k3)),...
			 numG,YprimeEF);
		vq_atk4 = RK_vq(delta(n,:)+dt.*delta_k3,Yg,RK_EGDQ(numG,delta(n,:)+dt.*delta_k3,...
			 RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt.*edd_k3,ed(n,:)+dt.*ed_k3),...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt.*eqq_k3,eq(n,:)+dt.*eq_k3)),...
			 numG,YprimeEF);
		delta_k4 = RK_delta(w(n,:)+dt.*w_k3,dt,f,n);
		w_k4 = RK_w(w(n,:)+dt.*w_k3,...
			 dt,n,f,H,...
			 Pm(n,:)+dt.*Pm_k3,...
			 RK_Pe(RK_egd(Kq,xqq,xqqq,xl,edd(n,:)+dt.*edd_k3,ed(n,:)+dt.*ed_k3),...
			 id_atk4,...
			 RK_egq(Kd,xdd,xddd,xl,eqq(n,:)+dt.*eqq_k3,eq(n,:)+dt.*eq_k3),...
			 iq_atk4,w0,xddd,xqqq,w(n,:)+dt.*w_k3),D,numG);
		eq_k4 = RK_eq(numG,Td,ef(n,:)+dt.*ef_k3,xd,xdd,xddd,xl,Kd,eqq(n,:)+dt.*eqq_k3,...
			 eq(n,:)+dt.*eq_k3,w(n,:)+dt.*w_k3,id_atk4,n,w0);
		eqq_k4 =RK_eqq(numG,Tdd,Kd,eq(n,:)+dt.*eq_k3,eqq(n,:)+dt.*eqq_k3,w(n,:)+dt.*w_k3,xdd,xl,...
			 id_atk4,n,w0);
		ed_k4 = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd(n,:)+dt.*edd_k3,w(n,:)+dt.*w_k3,iq_atk4,...
			 n,ed(n,:)+dt.*ed_k3,w0);
		edd_k4 =RK_edd(numG,Tqq,Kq,edd(n,:)+dt.*edd_k3,ed(n,:)+dt.*ed_k3,w(n,:)+dt.*w_k3,xqq,xl,...
			 iq_atk4,n,w0);
		ef_k4 = RK_ef(TA,KA,vd_atk4,vq_atk4,ef0,V,n,ef(n,:)+dt.*ef_k3);
		Pm_k4 = RK_Pm(Pm(n,:)+dt.*Pm_k3,w(n,:)+dt.*w_k3,numG,TG,KG,w0,n,Pm0);
		%--------- k4 ---------
%}
		delta(n+1,:) = delta(n,:) + dt./6.*(delta_k(2,:)+2.*delta_k(3,:)+2.*delta_k(4,:)+delta_k(5,:));
		w(n+1,:) = w(n,:) + dt./6.*(w_k(2,:)+2.*w_k(3,:)+2.*w_k(4,:)+w_k(5,:));
		eq(n+1,:) = eq(n,:) + dt./6.*(eq_k(2,:)+2.*eq_k(3,:)+2.*eq_k(4,:)+eq_k(5,:));
		eqq(n+1,:) = eqq(n,:) + dt./6.*(eqq_k(2,:)+2.*eqq_k(3,:)+2.*eqq_k(4,:)+eqq_k(5,:));
		ed(n+1,:) = ed(n,:) + dt./6.*(ed_k(2,:)+2.*ed_k(3,:)+2.*ed_k(4,:)+ed_k(5,:));
		edd(n+1,:) = edd(n,:) + dt./6.*(edd_k(2,:)+2.*edd_k(3,:)+2.*edd_k(4,:)+edd_k(5,:));
		ed(n+1,:) = ed(n,:) + dt./6.*(ed_k(2,:)+2.*ed_k(3,:)+2.*ed_k(4,:)+ed_k(5,:));
		ef(n+1,:) = ef(n,:) + dt./6.*(ef_k(2,:)+2.*ef_k(3,:)+2.*ef_k(4,:)+ef_k(5,:));
		Pm(n+1,:) = Pm(n,:) + dt./6.*(Pm_k(2,:)+2.*Pm_k(3,:)+2.*Pm_k(4,:)+Pm_k(5,:));

		egd(n+1,:) = RK_egd(Kq,xqq,xqqq,xl,edd(n+1,:),ed(n+1,:));
		egq(n+1,:) = RK_egq(Kd,xdd,xddd,xl,eqq(n+1,:),eq(n+1,:));
		id(n+1,:) = RK_id(delta(n+1,:),Yg,RK_EGDQ(numG,delta(n+1,:),egd(n+1,:),egq(n+1,:)),numG);
		iq(n+1,:) = RK_iq(delta(n+1,:),Yg,RK_EGDQ(numG,delta(n+1,:),egd(n+1,:),egq(n+1,:)),numG);
		vd(n+1,:) = RK_vd(delta(n+1,:),Yg,RK_EGDQ(numG,delta(n+1,:),egd(n+1,:),egq(n+1,:)),numG,YprimeEF);
		vq(n+1,:) = RK_vq(delta(n+1,:),Yg,RK_EGDQ(numG,delta(n+1,:),egd(n+1,:),egq(n+1,:)),numG,YprimeEF);
		Pe(n+1,:) = RK_Pe(egd(n+1,:),id(n+1,:),egq(n+1,:),iq(n+1,:),w0,xddd,xqqq,w(n+1,:));


		delta_for_plot(n,1) = dt*(n-1);
		delta_for_plot(n,2) = delta(n,2) - delta(n,1);
	end
	%//////////////// for loop ////////////////////
	%/////////////////////////////////////////////
	plot(delta_for_plot(:,1),delta_for_plot(:,2)/pi*180)
	%csvwrite('delta2_delta1.csv',delta_for_plot);
end


%////////////////////////// differential equation //////////////////////////
%--------- (8.2) ---------
function [ddeltadt] = RK_delta(w,dt,f,n)
	w0 = 2*pi*f;
	ddeltadt = w - w0;
end
%--------- (8.2) ---------

%--------- (8.3) ---------
function [dwdt] = RK_w(w,dt,n,f,H,Pm,Pe,D,numG)
	dwdt = zeros(1,numG);
	w0 = 2*pi*f;
	dwdt = (w0./2./H.*(w0./w.*Pm - w0./w.*Pe - ...
		 D.*(w./w0 - 1))) ;
end
%--------- (8.3) ---------

%--------- (8.4) ---------
function [deqdt] = RK_eq(numG,Td,ef,xd,xdd,xddd,xl,Kd,eqq,eq,w,id,n,w0)
	deqdt = zeros(1,numG);
	deqdt = 1./Td.*(ef + (xd - xdd).*(xdd - xddd)./((xdd - ...
		 xl).^2).*Kd.*eqq - ...
		 (1+(xd-xdd).*(xdd-xddd)./((xdd-xl).^2)).*eq - ...
		 w./w0.*(xd-xdd).*(xddd-xl)./(xdd-xl).*id);
end
%--------- (8.4) ---------

%--------- (8.5) ---------
function [deqqdt] = RK_eqq(numG,Tdd,Kd,eq,eqq,w,xdd,xl,id,n,w0)
	deqqdt = zeros(1,numG);
	deqqdt = ...
		 -1./Tdd./Kd.*(Kd.*eqq-eq+w./w0.*(xdd-xl).*id);
end
%--------- (8.5) ---------

%--------- (8.6) ---------
function [deddt] = RK_ed(numG,Tq,xq,xqq,xqqq,xl,Kq,edd,w,iq,n,ed,w0)
	deddt = zeros(1,numG);
	deddt = -1./Tq.*(-(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2).*Kq.*edd...
		 +(1+(xq-xqq).*(xqq-xqqq)./((xqq-xl).^2)).*ed ...
		 -w./w0.*(xq-xqq).*(xqqq-xl)./(xqq-xl).*iq);
end
%--------- (8.6) ---------
%--------- (8.7) ---------
function [dedddt] = RK_edd(numG,Tqq,Kq,edd,ed,w,xqq,xl,iq,n,w0)
	dedddt = zeros(1,numG);
	dedddt = -1./Tqq./Kq.*(Kq.*edd-ed-w./w0.*(xqq-xl).*iq);
end
%--------- (8.7) ---------

%--------- (8.115) ---------
function [defdt] = RK_ef(TA,KA,vd,vq,ef0,V,n,ef)
	defdt = -1./TA.*(ef+KA.*sqrt(vd.^2+vq.^2)...
		 -(ef0+KA.*V(1,:))); %V hanyousei
end
%--------- (8.115) ---------

%--------- (8.116) ---------
function [dPmdt] = RK_Pm(Pm,w,numG,TG,KG,w0,n,Pm0)
	dPmdt = zeros(1,numG);
	dPmdt = -1./TG .* (Pm + KG ./ w0 .* w - (Pm0(1,:) + KG));
end
%--------- (8.116) ---------
%////////////////////////// differential equation //////////////////////////


function [id] = RK_id(delta,Yg,EGDQ,numG)
	IGDQ = Yg*EGDQ;
	for k = 1:numG
		id(k) = real((IGDQ(2*k-1)+j*IGDQ(2*k))*exp(j*(pi/2-delta(k))));
	end
end

function [iq] = RK_iq(delta,Yg,EGDQ,numG)
	IGDQ = Yg*EGDQ;
	for k = 1:numG
		iq(k) = imag((IGDQ(2*k-1)+j*IGDQ(2*k))*exp(j*(pi/2-delta(k))));
	end
end

function [vd] = RK_vd(delta,Yg,EGDQ,numG,YprimeEF)
	IGDQ = Yg*EGDQ;
	VGDQ = inv(YprimeEF)*IGDQ;
	for k = 1:numG
		vd(k) = real((VGDQ(2*k-1)+j*VGDQ(2*k))*exp(j*(pi/2-delta(k))));
	end
end

function [vq] = RK_vq(delta,Yg,EGDQ,numG,YprimeEF)
	IGDQ = Yg*EGDQ;
	VGDQ = inv(YprimeEF)*IGDQ;
	for k = 1:numG
		vq(k) = imag((VGDQ(2*k-1)+j*VGDQ(2*k))*exp(j*(pi/2-delta(k))));
	end
end

function [EGDQ] = RK_EGDQ(numG,delta,egd,egq)
	EGDQ = zeros(2*numG,1);
	EGD = zeros(numG,1);
	EGQ = zeros(numG,1);
	EGD = egd .* sin(delta) + egq .* cos(delta);
	EGQ = -egd .* cos(delta) + egq .* sin(delta);
	for k = 1:numG
		EGDQ(2*k-1) = EGD(k);
		EGDQ(2*k) = EGQ(k);
	end
end

function [egd] = RK_egd(Kq,xqq,xqqq,xl,edd,ed)
	egd = Kq.*(xqq-xqqq)./(xqq-xl).*edd + (xqqq-xl)./(xqq-xl).*ed;
end

function [egq] = RK_egq(Kd,xdd,xddd,xl,eqq,eq)
	egq = Kd.*(xdd-xddd)./(xdd-xl).*eqq + (xddd-xl)./(xdd-xl).*eq;
end

function [Pe] = RK_Pe(egd,id,egq,iq,w0,xddd,xqqq,w)
	Pe = (egd.*id+egq.*iq-(xddd-xqqq).*id.*iq).*w./w0;
end


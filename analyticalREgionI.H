// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void analyticalRegion::calculateMassFlux(){

	R0square=sqr(R0);

	K = 2*k_g/rho_l/c_pg*log(B+1);
	Info<<"The liquid density is "<<rho_l<<endl;
	K_ = K;
	B_ = B;
	Info<<"K is "<<K<<endl;
	Rsquare = R0square-K*Time_;
	R = sqrt(Rsquare);
	R_ = R;

	q_dot_pp_R0 = heatFlux_;
	q_dot_pp_R = q_dot_pp_R0*R0square/Rsquare;

	//m_dot_pp_R = B*q_dot_pp_R/(Y_O_infty*stoi_ratio*Delta_H+c_pg*(T_infty-T_boil));
	m_dot_pp_R = q_dot_pp_R/(L+c_pl*(T_boil-T_infty));

	m_dot_pp_R0 = m_dot_pp_R*Rsquare/R0square;

	//B = exp(m_dot_pp_R*R*c_pg/k_g)-1;



	Info<<"mass flux of the interface is  "<<m_dot_pp_R0<<" "<<heatFlux_<<endl;

	//m_dot_pp_R0 = k_g/c_pg/R0*log(B+1);

	massFlux_ = m_dot_pp_R0;

	if(R<0){
		Info<<"The burn-out time is "<<Time_<<endl;
		FatalError << "The droplet is burnt out" << exit(FatalError);
	}

	Info<<"The radius of the droplet is "<<R<<endl;
	Info<<"The K is "<<K<<endl;

	FSR = log(B+1)/log(stoi_ratio*Y_O_infty+1);
	FSR_ = FSR;

	Info<<"The FSR is "<<FSR<<endl;

}

void analyticalRegion::calculateTemperature(){
	scalar Part_1;
	//scalar Part_6;

	Part_1 = pow((B+1),(1-R/R0));
	//Part_1 = exp(-m_dot_pp_R*Rsquare/k_g*c_pg/R0);
	//Part_6 = c_pg*(T_infty-T_boil)+Delta_H*stoi_ratio*Y_O_infty+L;
	interfaceTemperature_ = L*(Part_1-1)/c_pg+T_boil;;
	//interfaceTemperature_ = (Part_6*Part_1-L)/c_pg+T_boil;
	if( interfaceTemperature_ > 800 ){ interfaceTemperature_ = 800; }
}

void analyticalRegion::calculateFuel(){
	scalar Part_2;

	Part_2 = 1-exp(-m_dot_pp_R*Rsquare/k_g*c_pg/R0);
	interfaceFuel_ = (1+stoi_ratio*Y_O_infty)*Part_2-stoi_ratio*Y_O_infty;
	if(interfaceFuel_<0){interfaceFuel_ = 0.01;}

}

void analyticalRegion::calculateOxygen(){
	scalar Part_5;

	Part_5 = exp(-m_dot_pp_R*Rsquare/k_g*c_pg/R0)-1;

	interfaceOxygen_ = (1+stoi_ratio*Y_O_infty)*Part_5/stoi_ratio+Y_O_infty;
}

void analyticalRegion::calculateAssimilation_FSR()
{
	scalar ln_Bplus1;

	ln_Bplus1 = observation_*log(stoi_ratio*Y_O_infty+1);
	B = exp(ln_Bplus1)-1;

	B_ = B;
	Info<<"The B number is "<<B<<endl;
}

void analyticalRegion::calculateAssimilation_dropletSize()
{

	scalar targetPara;

	//inverse Part

	scalar Part_3;
	scalar Part_4;
	scalar observedK;

	observedK = (R0square-observation_)/Time_;
	//rho_l = 2*k_g*log(B+1)/observedK/c_pg;

	B = exp(observedK*rho_l*c_pg/2/k_g)-1;

	//if(B < 0){ B = 0.005; }

	if(rho_l < 0){ rho_l = 684; }

	Part_3 = (Y_O_infty*Delta_H*stoi_ratio+c_pg*(T_infty-T_boil))/B-Delta_H;
	Part_4 = -observation_*log(std::abs(Part_3)/c_pg/(T_boil-Y_O_infty));
	targetPara = Part_4/Time_;

	targetPara_ = targetPara;

	//forward Part
	//B = (Y_O_infty*Delta_H*stoi_ratio+c_pg*(T_infty-T_boil))/(c_pg*(T_boil-Y_O_infty)*exp(-targetPara/Rsquare*Time_));


	Info<<"Assimilation droplet size "<<observedK<<" "<<targetPara<<" "<<B<<endl;

}

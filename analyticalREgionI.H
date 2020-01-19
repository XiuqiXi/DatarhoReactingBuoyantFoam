// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void analyticalRegion::calculateMassFlux(){

	R0square=sqr(R0);

	K = 2 * k_g / rho_l / c_pg * Foam::log(B + 1);
	Info<<"The liquid density is "<<rho_l<<endl;
	K_ = K;
	B_ = B;
	Info<<"K is "<<K<<endl;
	Rsquare = R0square-coeff * K * Time_;
	R = sqrt(Rsquare);
	R_ = R;

	q_dot_pp_R0 = heatFlux_;
	q_dot_pp_R = q_dot_pp_R0*R0square/Rsquare;

	m_dot_pp_R = B * q_dot_pp_R / (Y_O_infty * stoi_ratio * Delta_H + c_pg * (T_infty - T_boil));
	//m_dot_pp_R = q_dot_pp_R/(L+c_pl*(T_boil-T_infty));

	m_dot_pp_R0 = m_dot_pp_R * Rsquare / R0square;

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

	FSR = log(B+1) / log(stoi_ratio * Y_O_infty + 1);
	FSR_ = FSR;

	Info<<"The FSR is "<<FSR<<endl;

}

void analyticalRegion::calculateStreamwiseVelocity(){
	if (similarVar>1){
		streamwiseVelocity_ = 0;
	}
	else{
		streamwiseVelocity_ = u_x*similarVar*pow((1-similarVar), 0.6666666);
	}
}

void analyticalRegion::calculateTransverseVelocity(){
	transverseVelocity_ = 0.002;
}

void analyticalRegion::calculateTemperature(){
	if (similarVar>1){
		interfaceTemperature_ = T_infty;
	}
	else{
		scalar J = pow((1-similarVar), 0.6666666);
		if (y<y_fl){
			interfaceTemperature_ = (T_boil-T_infty)*((1-D_c)*J+D_c)+T_infty;
		}
		else{
			interfaceTemperature_ = (T_boil-T_infty)*(D_c/r+1)*J+T_infty;
		}
	}
}

void analyticalRegion::calculateFuel(){
	if (similarVar>1){
		interfaceFuel_ = 0;
	}
	else{
		if (y<y_fl){
			scalar J = pow((1-similarVar), 0.6666666);
			interfaceFuel_ = Y_fw*((1+r)*J-r);
		}
		else{
			interfaceFuel_ = 0;
		}
	}
}

void analyticalRegion::calculateOxygen(){
	if (similarVar>1){
		interfaceOxygen_ = Y_O_infty;
	}
	else{
		scalar J = pow((1-similarVar), 0.6666666);
		if (y<y_fl){
			interfaceOxygen_ = 0;
		}
		else{
			interfaceOxygen_ = Y_O_infty*(1-(1+1/r)*J);
		}
	}
}

void analyticalRegion::calculateAssimilation_FSR()
{
	scalar ln_Bplus1;
	scalar B_temp;
	scalar observedK;

	ln_Bplus1 = observation_ * log(stoi_ratio * Y_O_infty + 1);
	B_temp = exp(ln_Bplus1) - 1;

	observedK = 2 * k_g / rho_l / c_pg * log(B_temp + 1);

	coeff = observedK / K;

	if(coeff<0){coeff = 1;}

	Info<<"Coefficient is "<<coeff<<endl;

}

void analyticalRegion::calculateAssimilation_dropletSize()
{

	scalar observedK;

	observedK = (R0square-observation_) / Time_;
	coeff = observedK / K;

	if(coeff<0){coeff = 1;}

	Info<<"Coefficient is "<<coeff<<endl;

}

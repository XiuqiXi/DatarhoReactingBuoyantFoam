// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void analyticalRegion::calculateStreamwiseVelocity(){
	if (similarVar>1){
		streamwiseVelocity_ = 0;
	}
	else{
		streamwiseVelocity_ = u_x*similarVar*pow((1-similarVar), 0.6666666);
	}
}

void analyticalRegion::calculateTransverseVelocity(){
	scalar v_sureface =  3*Foam::log(B+1)*k_g/2/c_pg/rho_boil/delta_1/pow(x, 0.25);
	if (similarVar>1){
		streamwiseVelocity_ = 0;
	}
	else{
		transverseVelocity_ = v_sureface+((u_1*(delta_1*pow(x, 0.25) - zeta)*\
		Foam::sqrt(1 - zeta/(delta_1*pow(x, 0.25))+VSMALL)*(12*delta_1*delta_1*Foam::sqrt(x)+18*delta_1*pow(x, 0.25)*zeta + 5*zeta*zeta))/\
		(140*delta_1*delta_1*x+VSMALL))-3*delta_1*u_1/35/pow(x, 0.25);
;
	}
	if (transverseVelocity_<0){transverseVelocity_ = 0;}
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

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

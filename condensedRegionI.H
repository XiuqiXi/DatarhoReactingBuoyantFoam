// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void condensedRegion::pyrolysisLengthIncrease_thermally_thick(){
  pyrolysisLengthIncrease_ = timeInterval_*preheatLength_/(14*constant::mathematical::pi/4*k*rho*c*pow((T_ignition-T_infty)/heatFlux_,2));
}

void condensedRegion::pyrolysisLengthIncrease_thermally_thin(){
  pyrolysisLengthIncrease_ = timeInterval_*preheatLength_/(rho*c*d*(T_ignition-T_infty)/heatFlux_);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void condensedRegion::pyrolysisLengthIncrease_thermally_thick(){
  pyrolysisLengthIncrease_ = preheatLength_/(4*constant::mathematical::pi*k*rho*c/pow((T_ignition-T_infty)/heatFlux_,2));
}

void condensedRegion::pyrolysisLengthIncrease_thermally_thin(){
  pyrolysisLengthIncrease_ = preheatLength_/(rho*c*d*(T_ignition-T_infty)/heatFlux_);
}

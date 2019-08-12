		// find ID of patch
		label patchID = mesh.boundaryMesh().findPatchID(interfaceName);
		// check patch has been found
		if(patchID == -1)
		{
				FatalError << "patch not found!" << exit(FatalError);
		}


		//calculate the heat flux of the inerface
		volScalarField& he = thermo.he();
		const volScalarField& alpha = thermo.alpha();

		Info<<"The maximun of alpha is "<<gMax(alpha.boundaryField()[patchID])<<"The minimun of alpha is "<<gMin(alpha.boundaryField()[patchID])<<endl;

		//the heat flux of each surfaces of each volumes
		heatFlux = fvc::interpolate(alpha)*fvc::snGrad(he);

		//change the surface heat flux as the volume(vol) heat flux
		volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

		const surfaceScalarField::Boundary& heatFluxBf = heatFlux.boundaryField();

		//calculate the mass flux of the interface

		scalar averageHeatBC = 0;
		scalar averageMassBC = 0;

		scalar minHfp = 0;
		scalar maxHfp = 0;
		scalar meanHfp = 0;

		//find the cylinder patch to calculate the heat flux
		wallHeatFluxBf[patchID] = heatFluxBf[patchID];
		wallHeatFlux = -wallHeatFlux;
		const scalarField& hfp = wallHeatFlux.boundaryField()[patchID];
		//the average heat
		minHfp = gMin(hfp);
		maxHfp = gMax(hfp);
		meanHfp = (minHfp+maxHfp)/2;
		averageHeatBC = meanHfp;

		Info<<endl;
		Info<<"The minimun value of "<<patchID<< " the heat flux is "<<minHfp<<endl;
		Info<<"The maximun value of "<<patchID<<" the heat flux is "<<maxHfp<<endl;
		Info<<endl;


		Info<<"The average value of the heat flux is "<<averageHeatBC<<endl;
		Info<<endl;

		/*
		Calculate the mass flux(slove the diffusion problem)
		*/

		//claim yhe Droplet class


		analyticalRegion Droplet(averageHeatBC,mesh.time().value(),originalDropletSize.value(),assimilationObservation);


		Droplet.calculateMassFlux();
		averageMassBC = Droplet.massFlux_;
		Info<<"The original size of the droplet is "<<Droplet.R0<<endl;
		/*
		end of calculating the mass flux
		*/

		Info<<"The average mass flux is "<<averageMassBC<<endl;

		//The injection velocity of the boundary
		scalar velocity_inject=averageMassBC/1;

		//mapping the calculated value on the velocity boundary condition

		Info<<endl;

		//Cast the boundary patch into the surfaceNormalFixed object
		surfaceNormalFixedValueFvPatchVectorField& Upatch = refCast<surfaceNormalFixedValueFvPatchVectorField>(U.boundaryFieldRef()[patchID]);

		//update the normal velocity
		//if(velocity_inject < 0){velocity_inject = 0.001;}
		if(mesh.time().value()<0.05){velocity_inject = 0.00;}
		Upatch.refValue_=-velocity_inject;

		//update the temperature

		scalar interfaceTemperature;
		Droplet.calculateTemperature();
		interfaceTemperature = Droplet.interfaceTemperature_;
		if(mesh.time().value()<0.05){interfaceTemperature = 375;}
		Info<<"interfaceTemperature is "<<interfaceTemperature<<endl;

		T.boundaryFieldRef()[patchID] == interfaceTemperature;

		scalar interfaceFuel;
		Droplet.calculateFuel();
		interfaceFuel = Droplet.interfaceFuel_;
		if(mesh.time().value()<0.1){interfaceFuel = 0.2;}
		Info<<"interfaceFuel is "<<interfaceFuel<<endl;
		composition.Y("C7H16").boundaryFieldRef()[patchID] == interfaceFuel;

		Info<<"B number is "<<analyticalRegion::B<<endl;

		//Info<<"C3H8 boundary condition is "<<fuel<<endl;

		//write the mass flux data into AverageMassFlux file
		std::ofstream file_velocity;

		file_velocity.open("InjectionVelocity", std::ofstream::out | std::ofstream::app);

		file_velocity<<runTime.timeName()<<token::TAB<<velocity_inject<<std::endl;

		file_velocity.close();

		//write the droplet properties
		std::ofstream file_droplet;

		file_droplet.open("DropletProperties",std::ofstream::out | std::ofstream::app);

		file_droplet<<runTime.timeName()<<token::TAB<<Droplet.R_<<token::TAB<<interfaceTemperature<<token::TAB<<interfaceFuel<<token::TAB<<Droplet.B_<<token::TAB<<Droplet.K_<<std::endl;

		file_droplet.close();

		//write the heat flux
		std::ofstream file_heatFlux;

		file_heatFlux.open("HeatFlux", std::ofstream::out | std::ofstream::app);

		file_heatFlux<<runTime.timeName()<<token::TAB<<maxHfp<<token::TAB<<minHfp<<token::TAB<<averageHeatBC<<std::endl;

		file_heatFlux.close();

		//write the mass flux
		std::ofstream file_massFlux;

		file_massFlux.open("MassFlux", std::ofstream::out | std::ofstream::app);

		file_massFlux<<runTime.timeName()<<token::TAB<<averageMassBC<<token::TAB<<Droplet.B_<<std::endl;

		file_massFlux.close();

		//Assimilation Loop

		if((std::abs(mesh.time().value() - data[assimilationTimeIndex][0])<0.001) && assimilationActivation == true && assimilationDataType == "dropletSize")
		{

			assimilationObservation = data[assimilationTimeIndex][1];

			Droplet.calculateAssimilation_dropletSize();

			std::ofstream file;

			file.open("assimilationData", std::ofstream::out | std::ofstream::app);

			file<<runTime.timeName()<<token::TAB<<assimilationTimeIndex<<token::TAB<<data[assimilationTimeIndex][0]<<token::TAB<<data[assimilationTimeIndex][1]<<token::TAB<<std::endl;

			file.close();

			assimilationTimeIndex++;
		}

		if((std::abs(mesh.time().value() - data[assimilationTimeIndex][0])<0.001) && assimilationActivation == true && assimilationDataType == "FSR")
		{

			assimilationObservation = data[assimilationTimeIndex][1];

			Droplet.calculateAssimilation_FSR();

			std::ofstream file;

			file.open("assimilationData", std::ofstream::out | std::ofstream::app);

			file<<runTime.timeName()<<token::TAB<<assimilationTimeIndex<<token::TAB<<data[assimilationTimeIndex][0]<<token::TAB<<data[assimilationTimeIndex][1]<<std::endl;

			file.close();

			assimilationTimeIndex++;
		}


		Info<<"The index is "<<	assimilationTimeIndex <<endl;
		Info<<"The time is "<<	data[assimilationTimeIndex][0] <<endl;
		Info<<"The observation is "<<data[assimilationTimeIndex][1] <<endl;

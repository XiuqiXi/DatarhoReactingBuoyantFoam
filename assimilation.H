		// find ID of patch
		label patchID = mesh.boundaryMesh().findPatchID(interfaceName);
		// check patch has been found
		if(patchID == -1)
		{
				FatalError << "patch not found!" << exit(FatalError);
		}

		scalar pyrolysisLowerPosition = pyrolysisStartPosition.value();
		scalar pyrolysisUpperPosition = pyrolysisLowerPosition+originalPyrolsisLength.value();


		//calculate the heat flux of the inerface
		volScalarField& he = thermo.he();
		const volScalarField& alpha = thermo.alpha();

		//the heat flux of each surfaces of each volumes
		heatFlux = fvc::interpolate(alpha)*fvc::snGrad(he);

		//change the surface heat flux as the volume(vol) heat flux
		volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

		const surfaceScalarField::Boundary& heatFluxBf = heatFlux.boundaryField();

		//find the interface patch to calculate the heat flux
		wallHeatFluxBf[patchID] = heatFluxBf[patchID];
		wallHeatFlux = -wallHeatFlux;
		const scalarField& hfp = wallHeatFlux.boundaryField()[patchID];

		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			//Find the location of preheat ergion
			/*if (hfp[faceI]>10000){
				Info<<"The location of preheat region is "<<y<<endl;
				break;
			}*/

			Info<<hfp[faceI]<<endl;
			//Info<<y<<endl;

		}


		//Define the velocity boundary condition object
		fvPatchVectorField& U_patch = refCast<fvPatchVectorField>(U.boundaryFieldRef()[patchID]);

		//Cast calculated velocity to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;
			//Info<<"Test"<<endl;

			if (y<pyrolysisLowerPosition){
				U_patch[faceI] = vector(0, 0, 0);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Info<<"value is "<<nearWallRegion.testValue<<endl;
				//Get the velocity
				nearWallRegion.calculateStreamwiseVelocity();
				scalar streamwiseVelocity = nearWallRegion.streamwiseVelocity_;
				nearWallRegion.calculateTransverseVelocity();
				scalar transverseVelocity = nearWallRegion.transverseVelocity_;
				Info<<"x is "<<x_original<<endl;
				Info<<"y is "<<y_original<<endl;
				Info<<"transverseVelocity is "<<transverseVelocity<<endl;
				Info<<"streamwiseVelocity is "<<streamwiseVelocity<<endl;
				//Set the velocity
				U_patch[faceI] = ( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition)\
				 ? vector(transverseVelocity, streamwiseVelocity, 0) : vector(0, 0, 0));
			}

		}
		Info<<"The velocity changes"<<endl;

		//update the temperature
		//Define the temperature boundary condition object
		fvPatchScalarField& T_patch = refCast<fvPatchScalarField>(T.boundaryFieldRef()[patchID]);

		//Cast calculated temperature to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;

			if (y<pyrolysisLowerPosition){
				T_patch[faceI] = scalar(293);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Get the temperature
				nearWallRegion.calculateTemperature();
				scalar interfaceTemperature = nearWallRegion.interfaceTemperature_;
				//Set the temperature
				T_patch[faceI] = scalar( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition) ? interfaceTemperature : 293 );
			}

		}
		Info<<"The temperature changes"<<endl;

		//update the fuel
		//Define the fuel boundary condition object
		fvPatchScalarField& fuel_patch = refCast<fvPatchScalarField>(composition.Y("C7H16").boundaryFieldRef()[patchID]);

		//Cast calculated fuel to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;

			if (y<pyrolysisLowerPosition){
				fuel_patch[faceI] = scalar(0);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Get the fuel
				nearWallRegion.calculateFuel();
				scalar interfaceFuel = nearWallRegion.interfaceFuel_;
				//Set fuel
				fuel_patch[faceI] = scalar( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition) ? interfaceFuel : 0 );
			}
		}
		Info<<"The fuel changes"<<endl;


		//update the oxygen
		//Define the oxygen boundary condition object
		fvPatchScalarField& oxygen_patch = refCast<fvPatchScalarField>(composition.Y("O2").boundaryFieldRef()[patchID]);

		//Cast calculated oxygen to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y-pyrolysisLowerPosition;

			if (y<pyrolysisLowerPosition){
				oxygen_patch[faceI] = scalar(0.23);
			}
			else{
				//Construct the analytical region
				analyticalRegion nearWallRegion(x_original, y_original);
				//Get the oxygen
				nearWallRegion.calculateOxygen();
				scalar interfaceOxygen = nearWallRegion.interfaceOxygen_;
				//Set oxygen
				oxygen_patch[faceI] = scalar( (y<pyrolysisUpperPosition)&&(y>pyrolysisLowerPosition) ? interfaceOxygen : 0.23 );
			}
		}
		Info<<"The oxygen changes"<<endl;



		Info<<"B number is "<<analyticalRegion::B<<endl;




		Info<<"The index is "<<	assimilationTimeIndex <<endl;
		Info<<"The time is "<<	data[assimilationTimeIndex][0] <<endl;
		Info<<"The observation is "<<data[assimilationTimeIndex][1] <<endl;

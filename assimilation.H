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


		Info<<"The average value of the heat flux  is "<<averageHeatBC<<endl;
		Info<<endl;



		//Define the velocity boundary condition object
		fvPatchVectorField& U_patch = refCast<fvPatchVectorField>(U.boundaryFieldRef()[patchID]);

		//Cast calculated velocity to each cell on the patch
		forAll(mesh.boundaryMesh()[patchID].faceCentres(), faceI)
		{
			//Get the coordinate
			const scalar x = mesh.boundaryMesh()[patchID].faceCentres()[faceI].x();
			const scalar y = mesh.boundaryMesh()[patchID].faceCentres()[faceI].y();
			scalar x_original = x+0.000002;
			scalar y_original = y+0.002;
			Info<<"Test"<<endl;

			//Construct the analytical region
			analyticalRegion nearWallRegion(x_original, y_original);
			Info<<"value is "<<nearWallRegion.testValue<<endl;
			//Get the velocity
			nearWallRegion.calculateStreamwiseVelocity();
			scalar streamwiseVelocity = nearWallRegion.streamwiseVelocity_;
			nearWallRegion.calculateTransverseVelocity();
			scalar transverseVelocity = nearWallRegion.transverseVelocity_;
			Info<<"x is "<<x_original<<endl;
			Info<<"y is "<<y_original<<endl;
			Info<<"V is "<<streamwiseVelocity<<endl;
			//Set the velocity
			U_patch[faceI] = ( (y<0.6) ? vector(streamwiseVelocity, transverseVelocity, 0) : vector(0, 0, 0));

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
			scalar y_original = y+0.002;

			//Construct the analytical region
			analyticalRegion nearWallRegion(x_original, y_original);
			//Get the temperature
			nearWallRegion.calculateTemperature();
			scalar interfaceTemperature = nearWallRegion.interfaceTemperature_;
			//Set the temperature
			T_patch[faceI] = scalar( (y<0.6) ? interfaceTemperature : 293 );

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
			scalar y_original = y+0.002;

			//Construct the analytical region
			analyticalRegion nearWallRegion(x_original, y_original);
			//Get the fuel
			nearWallRegion.calculateFuel();
			scalar interfaceFuel = nearWallRegion.interfaceFuel_;
			//Set fuel
			fuel_patch[faceI] = scalar( (y<0.6) ? interfaceFuel : 0 );

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
			scalar y_original = y+0.002;

			//Construct the analytical region
			analyticalRegion nearWallRegion(x_original, y_original);
			//Get the oxygen
			nearWallRegion.calculateOxygen();
			scalar interfaceOxygen = nearWallRegion.interfaceOxygen_;
			//Set oxygen
			oxygen_patch[faceI] = scalar( (y<0.6) ? interfaceOxygen : 0.23 );

		}
		Info<<"The oxygen changes"<<endl;



		Info<<"B number is "<<analyticalRegion::B<<endl;




		Info<<"The index is "<<	assimilationTimeIndex <<endl;
		Info<<"The time is "<<	data[assimilationTimeIndex][0] <<endl;
		Info<<"The observation is "<<data[assimilationTimeIndex][1] <<endl;

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
	Foam::analyticalRegion


Description
	Analytical solutions to the fuel droplet combustion problems. The region provides the conditions for the interface.


SourceFiles
    analyticalRegion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_analyticalRegion_H
#define Foam_analyticalRegion_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class analyticalRegion Declaration
\*---------------------------------------------------------------------------*/

class analyticalRegion
{


    //Static Menmbers

    //static scalar B;

    // Private Member
    scalar heatFlux_;
    scalar Time_;
    scalar observation_;


    scalar R0square;
    scalar R;
	  scalar Rsquare;

	const scalar Delta_H = 44926000;
	const scalar stoi_ratio = 0.32;

	const scalar T_boil = 371.5;
	const scalar T_infty = 293.15;
	const scalar Y_O_infty = 0.233;

  scalar FSR;

	scalar K;

	const scalar k_g = 0.0926;
	static scalar rho_l;
	//scalar rho_g = 1;
  const scalar c_pl = 2240;
	const scalar c_pg = 4220;
	const scalar L = 316000;

	//p represents prime notation
	scalar q_dot_pp_R;
	scalar q_dot_pp_R0;
	scalar m_dot_pp_R;
	scalar m_dot_pp_R0;


public:

   scalar B_;
   static scalar B;
	 scalar massFlux_;
	 scalar interfaceTemperature_;
	 scalar R_;
	 scalar interfaceFuel_;
   scalar interfaceOxygen_;
	 scalar R0;
   scalar targetPara_;
   scalar K_;
   scalar FSR_;

    // Constructors

        //- Construct from Time
        analyticalRegion(scalar heatFlux, scalar Time, scalar R0Origin, scalar observation);


	// Member functions

	void calculateMassFlux();
	void calculateTemperature();
	void calculateFuel();
  void calculateOxygen();
	void calculateAssimilation_FSR();
  void calculateAssimilation_dropletSize();

    //- Destructor
    ~analyticalRegion();
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "analyticalREgionI.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalRegion::analyticalRegion
(
    scalar heatFlux,
    scalar Time,
    scalar R0Origin,
    scalar observation
)
{
    heatFlux_ = heatFlux;
    Time_ = Time;
    R0 = R0Origin;
    observation_ = observation;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

analyticalRegion::~analyticalRegion()
{}


// ************************************************************************* //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// ************************************************************************* //

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
	Foam::condensedRegion


Description
	The Quintiere's model for condensed phase of the flame spread


SourceFiles
    condensedRegion.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_condensedRegion_H
#define Foam_condensedRegion_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class condensedRegion Declaration
\*---------------------------------------------------------------------------*/

class condensedRegion
{


    // Private Member
    const scalar k = 0.16;
    const scalar rho = 1200;
    const scalar c = 2100;
    scalar d = 0.001;

    const scalar T_ignition = 393.15;
    const scalar T_infty = 293;

    scalar timeInterval_;


public:

    scalar heatFlux_;
    scalar preheatLength_;
    scalar pyrolysisLengthIncrease_;


    // Constructors

    //- Construct from timeInterval
    condensedRegion(scalar heatFlux, scalar timeInterval, scalar preheatLength);


	   // Member functions

     //scalar calculateJ();
     void pyrolysisLengthIncrease_thermally_thick();
     void pyrolysisLengthIncrease_thermally_thin();

     //- Destructor
     ~condensedRegion();
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "condensedRegionI.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

condensedRegion::condensedRegion
(
    scalar heatFlux,
    scalar timeInterval,
    scalar preheatLength
)
{
    heatFlux_ = heatFlux;
    timeInterval_ = timeInterval;
    preheatLength_ = preheatLength;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

condensedRegion::~condensedRegion()
{}


// ************************************************************************* //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// ************************************************************************* //

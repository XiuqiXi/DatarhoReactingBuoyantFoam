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
	Analytical solutions to the boundary layer flame spread problems. The region provides the conditions for the interface.


SourceFiles
    analyticalRegion.H

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


    // Private Member

    scalar x;
    scalar y;

    const scalar Delta_H = 44926000;
    const scalar stoi_ratio = 0.32;

    const scalar T_boil = 371.5;
    const scalar T_infty = 293.15;
    const scalar Y_O_infty = 0.233;
    const scalar g = 9.81;

    const scalar rho_infty = 1.225;
    const scalar rho_boil = rho_infty*T_infty/T_boil;

    const scalar k_g = 0.0926;
    const scalar mu = 0.0000176;
    const scalar c_pl = 2240;
    const scalar c_pg = 4220;
    const scalar L = 316000;
    const scalar h_w = c_pg*(T_boil-T_infty);
    const scalar tau = h_w/L;

    scalar Y_fw = (B-stoi_ratio*Y_O_infty)/(1+B);
    scalar D_c = 1.2;
    scalar r = Y_O_infty*stoi_ratio/Y_fw;
    scalar J_fl = 1/(1+1/r);
    scalar T_fl = (T_boil-T_infty)*((1-D_c)*J_fl+D_c)+T_infty;
    scalar rho_fl = rho_infty*T_infty/T_fl;

    //calcuate assisting parameters
    scalar inter = 4/5*B*(Foam::pow(J_fl, 1.666666)-1)-2*(Foam::pow(J_fl, 0.666666)-1)*(B+tau);
    scalar M_1 = 48*inter*g*L/c_pg/T_infty;
    scalar M_2 = -48*mu*Foam::log(1+B)/B/rho_boil;
    scalar M_3 = 40*(B+1)*Foam::log(1+B)*k_g/B/c_pg/rho_infty;

    scalar delta_1 = Foam::pow((M_3*M_3-M_2*M_3)/M_1,0.25);

    scalar u_1 = M_3/Foam::pow(((M_3*M_3-M_2*M_3)/M_1), 0.5);

    scalar zeta_fl = delta_1*Foam::pow(x, 0.25)*(1-Foam::pow(J_fl, 0.666666));
    scalar zeta_delta = delta_1*Foam::pow(x, 0.25);

    scalar y_fl = 2*rho_boil*zeta_fl/(rho_fl+rho_boil+VSMALL);

    scalar Lambda_1 = (rho_fl-rho_boil)/rho_boil;
    scalar Lambda_2 = y_fl+VSMALL;
    scalar Lambda = Lambda_1/Lambda_2;


    scalar Gamma = Lambda*y_fl*y_fl/2+y_fl;
    scalar y_delta = (-2*Gamma*rho_boil+rho_infty*y_fl+rho_fl*y_fl+2*rho_boil*zeta_delta)/(rho_infty+rho_fl+VSMALL);

    scalar Pi_1 = (rho_infty-rho_fl)/(y_delta-y_fl+VSMALL)/rho_boil;
    scalar Pi_2 = (rho_infty-(rho_infty-rho_fl)/(y_delta-y_fl+VSMALL)*y_delta)/rho_boil;
    scalar Pi_3 = Gamma-(Pi_1*y_fl*y_fl/2+Pi_2*y_fl);

    scalar zeta = ((y<y_fl) ? (Lambda*y*y/2+y) : (Pi_1*y*y/2+Pi_2*y+Pi_3) );

    scalar u_x = u_1*Foam::pow(x, 0.5);
    scalar delta_x = delta_1*Foam::pow(x, 0.25);
    scalar similarVar = zeta/(delta_x+VSMALL);


public:

    scalar B_;
    static scalar B;
    scalar interfaceTemperature_;
    scalar interfaceFuel_;
    scalar interfaceOxygen_;
    scalar streamwiseVelocity_;
    scalar transverseVelocity_;

    scalar testValue = zeta;

    // Constructors

    //- Construct from Time
    analyticalRegion(scalar x_original, scalar y_original);


    // Member functions

    //scalar calculateJ();
    void calculateStreamwiseVelocity();
    void calculateTransverseVelocity();
  	void calculateTemperature();
  	void calculateFuel();
    void calculateOxygen();

    //- Destructor
    ~analyticalRegion();
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "analyticalREgionI.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalRegion::analyticalRegion
(
    scalar x_original,
    scalar y_original
)
{
    y = x_original;
    x = y_original;
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

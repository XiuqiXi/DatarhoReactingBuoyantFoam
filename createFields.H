#include "createRDeltaT.H"

Info<< "Creating reaction model\n" << endl;

autoPtr<combustionModels::rhoCombustionModel> reaction
(
    combustionModels::rhoCombustionModel::New(mesh)
);

rhoReactionThermo& thermo = reaction->thermo();
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

//pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the reaction model
reaction->setTurbulence(turbulence());

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

mesh.setFluxRequired(p_rgh.name());

#include "phrghEqn.H"

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);

surfaceScalarField heatFlux
(
    IOobject
    (
        "heatFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
	mesh,
    dimensionedScalar("heatFlux", dimMass/pow3(dimTime),scalar(0))
);

volScalarField wallHeatFlux
(
    IOobject
    (
        "wallHeatFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wallHeatFlux", dimMass/pow3(dimTime), 0.0)
);

volScalarField wallMassFlux
(
    IOobject
    (
        "wallMassFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wallMassFlux", dimMass/pow3(dimTime), 0.0)
);

IOdictionary analyticalRegionProperties
(
    IOobject
    (
        "analyticalRegionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar originalPyrolsisLength
(
	"originalPyrolsisLength",
	analyticalRegionProperties.lookup("originalPyrolsisLength")
);

dimensionedScalar pyrolysisStartPosition
(
	"pyrolysisStartPosition",
	analyticalRegionProperties.lookup("pyrolysisStartPosition")
);


bool assimilationActivation = readBool(analyticalRegionProperties.lookup("assimilationActivation"));


word interfaceName = word(analyticalRegionProperties.lookup("interface"));

analyticalRegion::B = readScalar(analyticalRegionProperties.lookup("initialBNumber"));

scalar pyrolysisLowerPosition = pyrolysisStartPosition.value();
scalar pyrolysisUpperPosition = pyrolysisLowerPosition+originalPyrolsisLength.value();

scalar preheatRegionUpperPosition;
scalar preheatLength;
scalar currentStepTime;
scalar timeInterval;
scalar previousStepTime = 0;

#include "createMRF.H"
#include "createRadiationModel.H"



Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.renameRho()
);









Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);






#include "compressibleCreatePhi.H"


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);






Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermoThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermoThermophysicalTransportModel::New(turbulence(), thermo)
);



dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"






Info << "Create field e_phi" << endl;

volScalarField e_phi
(
    IOobject
    (
        "e_phi",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);






Info<< "\nReading elProperties" << endl;

IOdictionary elProperties
(
    IOobject
    (
        "elProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar sigma0
(
    "sigma0",
    sqr(dimCurrent)*pow3(dimTime)/pow3(dimLength)/dimMass,
    elProperties.lookup("sigma0")
);

dimensionedScalar sigma1
(
    "sigma1",
    sqr(dimCurrent)*pow3(dimTime)/pow3(dimLength)/dimMass/dimTemperature,
    elProperties.lookup("sigma1")
);

dimensionedScalar sigma2
(
    "sigma2",
    sqr(dimCurrent)*pow3(dimTime)/pow3(dimLength)/dimMass/sqr(dimTemperature),
    elProperties.lookup("sigma2")
);

dimensionedScalar sigma3
(
    "sigma3",
    sqr(dimCurrent)*pow3(dimTime)/pow3(dimLength)/dimMass/pow3(dimTemperature),
    elProperties.lookup("sigma3")
);

dimensionedScalar sigma4
(
    "sigma4",
    sqr(dimCurrent)*pow3(dimTime)/pow3(dimLength)/dimMass/pow4(dimTemperature),
    elProperties.lookup("sigma4")
);




Info<< "Initialize buoyancy\n" << endl;



uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, 0)
);

uniformDimensionedScalarField pRef
(
    IOobject
    (
        "pRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimPressure, 0)
);

dimensionedScalar ghRef(-mag(g)*hRef);

volScalarField gh("gh", (g & mesh.C()) - ghRef);

surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

mesh.schemes().setFluxRequired(p_rgh.name());

Foam::pressureReference pressureReference(p,
        p_rgh,
        pimple.dict(),
        thermo.incompressible()
);

tmp<volVectorField> netForce;

netForce = new volVectorField
(
    IOobject
    (
        "netForce",
        runTime.name(),
        mesh
    ),
    fvc::reconstruct
    (
        (-ghf*fvc::snGrad(rho) - fvc::snGrad(p_rgh))
       *mesh.magSf()
    )
);






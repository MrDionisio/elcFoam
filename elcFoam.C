/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Application
    rhoPorousSimpleFoam

Description
    Steady-state solver for turbulent flow of compressible fluids, with
    implicit or explicit porosity treatment and optional sources.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermoThermophysicalTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "IOporosityModelList.H"

#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvcVolumeIntegrate.H"

#include "fvcSnGrad.H"
#include "fvcReconstruct.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "uniformDimensionedFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "out.H"
    #include "createZones.H"
    #include "initContinuityErrs.H"
    

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Adjust the time-step according to the solver maxDeltaT
//        adjustDeltaT(runTime, solver);

        runTime++;

        turbulence->predict();
        thermophysicalTransport->predict();

        // PIMPLE corrector loop
        
        while (pimple.loop())
        {
            #include "e_phiEqn.H"
            #include "UEqn.H"
            #include "EEqn.H"
            #include "pEqn.H"
        }
      
        turbulence->correct();
        thermophysicalTransport->correct();

        sigma = sigma4*pow4(thermo.T())+sigma3*pow3(thermo.T())+sigma2*sqr(thermo.T())+sigma1*thermo.T()+sigma0;
        C_p = thermo.Cp();
        dynamicViscosicy = thermo.mu();
        lambda_coef = thermo.kappa();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

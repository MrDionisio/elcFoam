{
    const volScalarField& psi = thermo.psi();
    rho = thermo.rho();
    rho.relax();

    fvVectorMatrix& UEqn = tUEqn.ref();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    const surfaceScalarField rhof(fvc::interpolate(rho));

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    tmp<volScalarField> rAtU
    (
        pimple.consistent()
      ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
      : tmp<volScalarField>(nullptr)
    );

    tmp<surfaceScalarField> rhorAtUf
    (
        pimple.consistent()
      ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
      : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        pimple.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        rhof*fvc::flux(HbyA)
//      + rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf)
    );

    MRF.makeRelative(rhof, phiHbyA);

    const bool adjustMass =
        mesh.schemes().steady() && adjustPhi(phiHbyA, U, p_rgh);

    const surfaceScalarField ghGradRhof(-ghf*fvc::snGrad(rho)*mesh.magSf());

    phiHbyA += rhorAUf*ghGradRhof;

    tmp<fvScalarMatrix> tp_rghEqn;

    {
        if (pimple.consistent())
        {
            const surfaceScalarField gradpf(fvc::snGrad(p_rgh)*mesh.magSf());
            phiHbyA += (rhorAAtUf - rhorAUf)*gradpf;
            HbyA += (rAAtU - rAU)*fvc::reconstruct(gradpf - ghGradRhof);
        }

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, rho, U, phiHbyA, rhorAAtUf, MRF);

//        fvc::makeRelative(phiHbyA, rho, U);

        fvScalarMatrix p_rghDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
          + fvc::div(phiHbyA)
//         ==
//            fvModels().source(psi, p_rgh, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            tp_rghEqn = p_rghDDtEqn - fvm::laplacian(rhorAAtUf, p_rgh);
            fvScalarMatrix& p_rghEqn = tp_rghEqn.ref();

            p_rghEqn.setReference
            (
                pressureReference.refCell(),
                pressureReference.refValue()
            );

            p_rghEqn.solve();
        }
    }

    const fvScalarMatrix& p_rghEqn = tp_rghEqn();

    phi = phiHbyA + p_rghEqn.flux();

    // Calculate and relax the net pressure-buoyancy force
    netForce.ref().relax
    (
        fvc::reconstruct((ghGradRhof + p_rghEqn.flux()/rhorAAtUf)),
        p_rgh.relaxationFactor()
    );

    // Correct the momentum source with the pressure gradient flux
    // calculated from the relaxed pressure
    U = HbyA + rAAtU*netForce();
    U.correctBoundaryConditions();

    p = p_rgh + rho*gh + pRef;

    // For steady compressible closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (adjustMass && !thermo.incompressible())
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh - pRef;
        p_rgh.correctBoundaryConditions();
    }

    // Optionally relax pressure for the thermophysics
    p.relax();

}

/*{
    const volScalarField& psi = thermo.psi();

    tmp<volVectorField> tHbyA;
    tHbyA = constrainHbyA(trAU()*UEqn.H(), U, p);
    volVectorField& HbyA = tHbyA.ref();

    tUEqn.clear();

    bool closedVolume = false;

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    );
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    closedVolume = adjustPhi(phiHbyA, U, p);

    while (pimple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;

        {
            tpEqn =
            (
                fvm::laplacian(rho*trAU(), p)
              + fvModels.source(psi, p, rho.name())
             ==
                fvc::div(phiHbyA)
            );
        }

        fvScalarMatrix& pEqn = tpEqn.ref();

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    {
        U = HbyA - trAU()*fvc::grad(p);
    }

    U.correctBoundaryConditions();

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();
    }

    rho = thermo.rho();
    rho.relax();
}
*/

    // Construct the Momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevTau(U)
     ==
        fvModels.source(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    tmp<volScalarField> trAU;

    if (pimple.momentumPredictor())
    {
            solve
            (
                UEqn
             ==
                netForce()
            );
    }

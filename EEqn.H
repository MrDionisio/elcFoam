{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(phi, he)
      + thermophysicalTransport->divq(he)
      - sigma * (fvc::grad(e_phi) & fvc::grad(e_phi)) 
     ==
        fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();
}

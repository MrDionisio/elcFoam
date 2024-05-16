{
    fvScalarMatrix e_phiEqn
    (
        fvm::laplacian(e_phi)
    );

    e_phiEqn.relax();

    e_phiEqn.solve();
}

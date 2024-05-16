        volScalarField C_p 
        (
            IOobject
            (
                "C_p",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo.Cp()
        );


        volScalarField dynamicViscosicy 
        (
            IOobject
            (
                "dynamicViscosicy",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo.mu()
        );

        volScalarField lambda_coef 
        (
            IOobject
            (
                "lambda_coef",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermo.kappa()
            
        );

        volScalarField sigma
        (
            IOobject
            (
                "sigma",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            sigma4*pow4(thermo.T())+sigma3*pow3(thermo.T())+sigma2*sqr(thermo.T())+sigma1*thermo.T()+sigma0   
        );

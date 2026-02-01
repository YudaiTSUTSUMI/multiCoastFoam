scalarList deltaTList(SWRegions.size(),0);

boolList limitFrontVelocitySW(SWRegions.size());
PtrList<dimensionedScalar> hDrySW(SWRegions.size());

PtrList<volScalarField> hSW(SWRegions.size());
PtrList<volScalarField> h0SW(SWRegions.size());
PtrList<volScalarField> hTotalSW(SWRegions.size());
PtrList<volVectorField> hUSW(SWRegions.size());
PtrList<volVectorField> USW(SWRegions.size());
PtrList<volScalarField> wetdrySW(SWRegions.size());
PtrList<volScalarField> NSW(SWRegions.size());

PtrList<SWRiemannSolver> riemannSolverSW(SWRegions.size());
PtrList<SWWaveSource> waveSourceSW(SWRegions.size());

forAll(SWRegions, i)
{
    Info<< "*** Region " << SWRegions[i].name() << nl << endl;

    IOdictionary SWProperties
    (
        IOobject
        (
            "SWProperties",
            runTime.constant(),
            SWRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    limitFrontVelocitySW[i] = SWProperties.lookupOrDefault("limitFrontVelocity", false);

    hDrySW.set
    (
        i,
        new dimensionedScalar
        (
            "hDry",
            dimLength,
            SWProperties.lookupOrDefault("hDry", SMALL)
        )        
    );

    Info<< "Reading field h\n" << endl;
    hSW.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "h",
                runTime.timeName(),
                SWRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            SWRegions[i]
        )
    );
    
    Info<< "Reading field h0\n" << endl;
    h0SW.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "h0",
                runTime.timeName(),
                SWRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            SWRegions[i]
        )
    );
    
    hTotalSW.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "hTotal",
                runTime.timeName(),
                SWRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            hSW[i]+h0SW[i]
        )
    );

    Info<< "Reading field hU\n" << endl;
    hUSW.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "hU",
                runTime.timeName(),
                SWRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            SWRegions[i]
        )
    );
    
    //Info<< "Reading field U\n" << endl;
    USW.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                SWRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            SWRegions[i],
     		dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
        )
    );
    
    //Info<< "Reading field wetdry\n" << endl;
	wetdrySW.set
	(
		i,
        new volScalarField
        (
            IOobject
            (
                "wetdry",
                runTime.timeName(),
                SWRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pos(hSW[i]-hDrySW[i])
        )
	);
    
    NSW.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "N",
                runTime.timeName(),
                SWRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            SWRegions[i],
            dimensionedScalar("N", dimensionSet(0,0,0,0,0,0,0),0.0)
        )
    );
        
    riemannSolverSW.set
    (
    	i,
    	new SWRiemannSolver
    	(
			hSW[i],
			h0SW[i],
			hUSW[i],
			NSW[i],
			wetdrySW[i],
			g.value()
    	)
    );
    
    waveSourceSW.set
    (
    	i,
    	new SWWaveSource
    	(
			runTime,
			hSW[i],
			hUSW[i]			
    	)
    );
}


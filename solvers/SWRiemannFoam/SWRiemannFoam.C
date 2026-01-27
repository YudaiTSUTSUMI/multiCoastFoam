/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    SWRiemannFoam

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "SWRiemannSolver.H"
#include "SWWaveSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"    
    #include "createMesh.H"    
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
    Info<< "\nStarting time loop\n" << endl;    

    while (runTime.run())
    {    
        waveSource.update();
        
        Info << "compute SWRiemannSolver flux" << endl;
        riemannSolver.computeFlux();

        if (adjustTimeStep)
        {        
            runTime.setDeltaT
            (
                min
                (
                    maxCo*riemannSolver.deltaT(),
                    maxDeltaT
                )
            );        
        }

	    runTime++;
        Info<< "\n deltaT = " <<  runTime.deltaTValue() << endl;
        Info<< "\n Time = " << runTime.timeName() << nl << endl;
                        
        Info << "compute SWRiemannSolver residue and update" << endl;
        riemannSolver.computeResidue();
        
        h = h - runTime.deltaT() * riemannSolver.hResidue();
        hU = hU - runTime.deltaT() * riemannSolver.hUResidue();
        
        wetdry = pos(h-hDry);
        h = wetdry*h;
        hU = wetdry*hU;   
        hTotal = h + h0;

        volScalarField h2 = h;
        h2.max(hDry);
        U = hU/h2;
        
        if(limitFrontVelocity)
        {
            volScalarField front = 0*wetdry;

            forAll(mesh.owner(), facei)
            {
                const label own = mesh.owner()[facei];
                const label nei = mesh.neighbour()[facei];

                if(wetdry[own] == 1 && wetdry[nei] == 0)
                {
                    front[own] = 1;
                }
                else if(wetdry[own] == 0 && wetdry[nei] == 1)
                {
                    front[nei] = 1;
                }
            }

            volScalarField UWet = wetdry*(1-front)*mag(U);
            scalar UMax = gMax(UWet);

            Info << "\n limited velocity: " << UMax << nl << endl; 

            forAll(U, celli)
            { 
                scalar magU = mag(U[celli]);
                if (front[celli] == 1.0 &&  magU > UMax)
                { 
                    U[celli] *= UMax / magU;
                    hU[celli] = h[celli]*U[celli];
                }
            }
        }        

        h.correctBoundaryConditions();
        hU.correctBoundaryConditions();
	    	    
        runTime.write();

        runTime.printExecutionTime(Info);
        Info << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

forAll(SWRegions, i)
{
	SWRiemannSolver& riemannSolver = riemannSolverSW[i];
	SWWaveSource& waveSource = waveSourceSW[i];  
	
	waveSource.update();
	
	Info << "compute " << SWNames[i] << " SWRiemannSolver flux" << endl;
    riemannSolver.computeFlux();
    
    Info << "compute " << SWNames[i] << " SWRiemannSolver residue and update" << endl;
    riemannSolver.computeResidue();
    
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
}    		        
  

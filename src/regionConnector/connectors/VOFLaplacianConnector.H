#ifndef VOFLaplacianConnector_H
#define VOFLaplacianConnector_H

#include "fvCFD.H"

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct VOFRegionPatchData
{
    const fvMesh* VOFMesh_;
    
    scalar patchID_;
    
    volScalarField* p_;
    
    labelList globalIDs_;
};

struct laplacianRegionData
{
    const fvMesh* laplacianMesh_;
    
    scalar patchID_;
    
    fixedValueFvPatchField<scalar>* p_;
    
    labelList globalIDs_;
};

class VOFLaplacianConnector
{    
    vector g_;
    
    word VOFPatchName_;
    word laplacianPatchName_;
            
    VOFRegionPatchData VOFData_;
    laplacianRegionData laplacianData_;
    
    labelList VOFIDsForLaplacian_; //global
    
    
    label findNearestFace
    (
    	const vector position,
		const vectorList& positionList
    );
    
    template<class templateList>
    templateList generateGlobalList
    (
        const templateList& list
    );
    
    labelList generateGlobalIDs
    (
        label localSize
    );
    
    
public:
    // Constructors
    VOFLaplacianConnector
    (
        const vector g,
        const dictionary& connectorDict,
        const fvMesh& VOFRegion,
        const fvMesh& laplacianRegion
    );

    //- Destructor
    ~VOFLaplacianConnector()
    {}
    
    void update();
    
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


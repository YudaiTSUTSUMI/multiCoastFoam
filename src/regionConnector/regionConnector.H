#ifndef regionConnector_H
#define regionConnector_H

#include "fvCFD.H"

#include "mixedFvPatchField.H"

#include "SWVOFFixedConnector.H"
#include "SWVOFArbitraryConnector.H"
#include "VOFLaplacianConnector.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionConnectors
{
    struct SWRegionData
    {
        const fvMesh* SWMesh_;
        
        volScalarField* cellType_; // 0, 1, 2, 3       
        surfaceScalarField* faceType_; // 0 or 1
        
        surfaceScalarField* hBound_;
        surfaceVectorField* hUBound_;
        
        PtrList<SWVOFFixedConnector> SWVOFFixedConnectorList_;
        PtrList<SWVOFArbitraryConnector> SWVOFArbitraryConnectorList_;
    };    
    
    // Private data
    vector g_;
    dictionary regionConnectProperties_;
    PtrList<SWVOFFixedConnector> SWVOFFixedConnectorList_;
    PtrList<SWVOFArbitraryConnector> SWVOFArbitraryConnectorList_;
    PtrList<VOFLaplacianConnector> VOFLaplacianConnectorList_;
    
    PtrList<SWRegionData> SWRegionData_;    
    
    // Private function
    void setSWRegionData
    (
        SWRegionData& SWData
    );
    
    void updateSWRegionData
    (
        SWRegionData& SWData
    );
            
public:

    // Constructors
    regionConnectors
    (
        const Time& runTime,
        const vector g,
        const PtrList<fvMesh>& SWRegions,
        const PtrList<fvMesh>& VOFRegions,
        const PtrList<fvMesh>& laplacianRegions,
        const wordList& SWRegionNames,
        const wordList& VOFRegionNames,
        const wordList& laplacianRegionNames
    );

    //- Destructor
    ~regionConnectors()
    {}
    
    void updateSWVOF();
    
    void updateVOFLaplacian();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


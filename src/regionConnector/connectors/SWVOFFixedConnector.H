#ifndef SWVOFFixedConnector_H
#define SWVOFFixedConnector_H

#include "fvCFD.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SWVOFFixedConnector
{
    // data 
    struct SWRegionData
    {
        //basic data
        const fvMesh* SWMesh_;
        
        volScalarField* h_;
        volVectorField* hU_;
        volScalarField* h0_;
        
        labelList cellType_;        
        labelList faceType_; // 0 or 1
        
        surfaceScalarField* hBound_;
        surfaceVectorField* hUBound_;
        
        //changed by mesh update
        labelList localFaceIDs_;
        
        labelList globalIDs_;
        
        List<labelList> VOFGlobalDonor_;
    };

    struct VOFRegionData
    {
        //basic data
        const fvMesh* VOFMesh_;
        
        volScalarField* alpha_;
        volScalarField* rho_;
        volVectorField* U_;
        
        scalar patchID_;
        word patchName_;
        
        mixedFvPatchField<scalar>* alphaBound_;
        fvPatchField<scalar>* rhoBound_;
        mixedFvPatchField<vector>* UBound_;
        fvsPatchField<scalar>* phiBound_;
        
        labelList globalIDs_; // local size, p0:[0, 1] p1:[2, 3, 4]
    };
    
    
    // private menber
    vector g_;
    
    SWRegionData SWRegionData_;
    
    VOFRegionData VOFRegionData_;
       
    vector VOFC_;
    
    // prvate function
    
        template<class templateList>
        templateList generateGlobalList
        (
            const templateList& list
        );
        
        labelList generateGlobalIDs
        (
            label localSize
        );
        
        template<class FieldType>
        List<FieldType> extractFieldSubset
        (
            const labelList& donorIDs,
            const UList<FieldType>& globalField
        );
            
    void updateConnection
    (
        scalar& h,
        scalar& hFraction,
       
        vector& hU,
        scalar& hUFraction,
        
        const scalar& h0,
        
        scalarField& alpha,
        scalarField& alphaFraction,
        
        vectorField& U,
        scalarField& UFraction,
        
        const vectorField& SfVOF,
        const vectorField& CVOF,
        const scalarField& deltaVOF,
        
        const scalar hVOF,
        const vector hUVOF,
        
        const scalar Fr
    );
    
    void storeBasicData
    (
        const dictionary& regionConnectorEntry,
        const fvMesh& SWRegion,
        const fvMesh& VOFRegion
    );
    
    void createConnection();
    
    void updateConnectBoundary();
    
public:
    // Constructors
    SWVOFFixedConnector
    (
        const vector g,
        const dictionary& regionConnectorEntry,
        const fvMesh& SWRegion,
        const fvMesh& VOFRegion
    );

    //- Destructor
    ~SWVOFFixedConnector()
    {}
    
    void update();
    
    const labelList& SWCellType() const
    {
        return SWRegionData_.cellType_;
    }
    
    const labelList& SWFaceType() const
    {
        return SWRegionData_.faceType_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


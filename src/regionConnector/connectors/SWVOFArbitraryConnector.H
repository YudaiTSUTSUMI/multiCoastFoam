#ifndef SWVOFArbitraryConnector_H
#define SWVOFArbitraryConnector_H

#include "fvCFD.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class SWVOFArbitraryConnector    
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
        labelList faceType_;
        
        surfaceScalarField* hBound_;
        surfaceVectorField* hUBound_;
        
        //changed by mesh update
        labelList localFaceIDs_;
        
        labelList globalIDs_;
        
        List<List<labelList>> VOFDonorProcCells_;
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
        
        vectorField SfBound_;
        vectorField CBound_;
        
        labelList globalIDs_; // local size, p0:[0, 1] p1:[2, 3, 4]
        
        //changed by mesh update
        labelList donorCell_;
        
        List<labelList> SWDonorProcCell_; //globalData
    };
    
    
    // private menber
    vector g_;
    
    SWRegionData SWRegionData_;
    
    VOFRegionData VOFRegionData_;
       
    vector VOFC_;
    
    bool updateConnection_;
    
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
    SWVOFArbitraryConnector
    (
        const vector g,
        const dictionary& regionConnectorEntry,
        const fvMesh& SWRegion,
        const fvMesh& VOFRegion
    );

    //- Destructor
    ~SWVOFArbitraryConnector()
    {}
    
    void update();
    
    bool updateConnection()
    {
        return updateConnection_;
    }
    
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


#ifndef SWWaveSource_H
#define SWWaveSource_H

#include "polyMesh.H"
#include "fvCFD.H"

#include "Time.H"
#include "csvArray.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class waveSourceData
{
    labelList hList_;
    List<scalarList> waveData_;
    scalar coeff_;    
    scalarList hInitList_;

    public:
    waveSourceData()
    :
        hList_(),
        waveData_()
    {}

    waveSourceData
    (
        const labelList& hList,
        const List<scalarList>& waveData,
        const scalar coeff,
        const scalarList hInitList
        
    )
    :
        hList_(hList),
        waveData_(waveData),
        coeff_(coeff),
        hInitList_(hInitList)
    {}
    
    inline labelList& hList() {return hList_;};
    inline List<scalarList>& waveData() {return waveData_;};
    inline scalar& coeff() {return coeff_;};
    inline scalarList& hInitList() {return hInitList_;};
};

class SWWaveSource
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;
        
        const Time& runTime_;
        
        volScalarField& h_;
        volVectorField& hU_;
        
        PtrList<waveSourceData> waveSources_;
       
public:

    // Constructors

        //- Construct from components
        SWWaveSource
        (
            const Time& runTime,
            volScalarField& h,
            volVectorField& hU
        );


    //- Destructor
    ~SWWaveSource()
    {}


    // Member Functions
    
    void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


#include "bulletBody.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::bulletBody::bulletBody
(
    const dictionary& bodyDict
)
:
    bodyDict_(bodyDict),
    dynamic_(false),
    P_(Zero),
    P0_(Zero),
    Q_(I),
    Q0_(I),
    initialP_(Zero),
    localCoM_(Zero),
    localCoM0_(Zero),
    v_(Zero),
    v0_(Zero),
    ev_(Zero),
    ev0_(Zero),
    a_(Zero),
    pi_(Zero),
    pi0_(Zero),
    tau_(Zero),
    rhoSolid_(0.0),
    massSolid_(0.0),
    inertiaSolid_(Zero),
    CoMSolid_(Zero),
    massFluid_(0.0),
    inertiaFluid_(Zero),
    inertiaTotal_(Zero),
    inertiaTotal0_(Zero),
    CoMFluid_(Zero)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar Foam::bulletBody::kinematicInfo::wrapToPi
(
    const scalar a
)
{
    scalar x = a;
    const scalar pi = constant::mathematical::pi;
    const scalar twoPi = 2.0*pi;

    while (x > pi)
    {
        x -= twoPi;
    }

    while (x < -pi)
    {
        x += twoPi;
    }

    return x;
}

scalar bulletBody::kinematicInfo::interpAngle
(
    const scalar a0,
    const scalar a1,
    const scalar theta
)
{
    const scalar da = wrapToPi(a1 - a0);
    return wrapToPi(a0 + theta*da);
}


void Foam::bulletBody::kinematicInfo::updateKinematicState
(
    scalar tNow
)
{
    const scalar degToRad = M_PI / 180.0;
    
    // ---- clamp after last sample ----
    if (tNow >= motionData_[motionData_.size()-1][0])
    {
        const label last = motionData_.size() - 1;
        linearState_ = vector(motionData_[last][1], motionData_[last][2], motionData_[last][3]);
        angularState_ = vector(motionData_[last][4], motionData_[last][5], motionData_[last][6]) * degToRad;
        return;
    }
    
    
    vector pos0;
    vector pos1;

    vector rotRad0;
    vector rotRad1;     
    
    scalar theta = 0.0;
    
    // ---- clamp before first sample ----
    if (tNow <= motionData_[0][0])
    {
        pos0 = linearState0_;        
        pos1 = vector(motionData_[0][1], motionData_[0][2], motionData_[0][3]);
        
        rotRad0 = angularState0_;
        rotRad1 = vector(motionData_[0][4], motionData_[0][5], motionData_[0][6]) * degToRad;
        
        const scalar t0 = 0.0;
        const scalar t1 = motionData_[0][0];
        
        if (t1 > t0 + SMALL)
        {
            theta = (tNow - t0)/(t1 - t0);
        }
        theta = max(min(theta, scalar(1)), scalar(0));
    }
    else
    {
        label k = 0;
        while (k < motionData_.size() - 2 && tNow >= motionData_[k+1][0])
        {
            k++;
        }
        
        pos0 = vector(motionData_[k][1],   motionData_[k][2],   motionData_[k][3]);
        pos1 = vector(motionData_[k+1][1], motionData_[k+1][2], motionData_[k+1][3]);

        rotRad0 = vector(motionData_[k][4],   motionData_[k][5],   motionData_[k][6]) * degToRad;
        rotRad1 = vector(motionData_[k+1][4], motionData_[k+1][5], motionData_[k+1][6]) * degToRad; 

        const scalar t0 = motionData_[k][0];
        const scalar t1 = motionData_[k+1][0];

        if (t1 > t0 + SMALL)
        {
            theta = (tNow - t0)/(t1 - t0);
        }
        theta = max(min(theta, scalar(1)), scalar(0));        
    }    

    linearState_ = (1.0 - theta)*pos0 + theta*pos1;

    angularState_ = vector
    (
        interpAngle(rotRad0.x(), rotRad1.x(), theta),  // roll
        interpAngle(rotRad0.y(), rotRad1.y(), theta),  // pitch
        interpAngle(rotRad0.z(), rotRad1.z(), theta)   // yaw
    );
}


// ************************************************************************* //

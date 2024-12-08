#include "cutFacePLICNew.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutFacePLICNew::cutFacePLICNew(const fvMesh& mesh)
:
    cutFace(mesh),
    mesh_(mesh),
    face_(),
    faceNormal_(Zero),
    eps_(1e-12),
    points_(mesh.points()),
    face_size_(0),
    pointStatus_(25),
    faceStatus_(-1)
{
}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
Foam::label Foam::cutFacePLICNew::map_status(const scalar &status1, const scalar &status2)
{
    if (status1 < -eps_)
        if (status2 < -eps_)
            return -1;
        else if (status2 > eps_)
            return 0;
        else return -2;
    else if (status1 > eps_)
        if (status2 > eps_)
            return 1;
        else if (status2 < -eps_)
            return 0;
        else return 2;
    else if (status2 < -eps_)
        return -2;
    else if (status2 > eps_)
        return 2;
    else return 3;
    
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //
void Foam::cutFacePLICNew::initCutFacePLICNew(
    const label cellI,
    const label faceI,
    const label owner,
    const scalarList& distanceList,
    const labelList& pLabels
)
{
    faceStatus_ = -2;
    
    // weird case that happens sometimes 
    if (owner != 0) {

    if (owner == cellI){
        order_ = false;}
    else{
        order_ = true;}

    } 
    else {
        if (((mesh_.faces()[faceI].centre(points_)-mesh_.C()[cellI]) & mesh_.faces()[faceI].areaNormal(points_) )> 0){ /////////// careful, if base not cell centre, change it to cell centre here
            order_ = false;}
        else {order_ = true;}
    }
    if (order_) face_ = mesh_.faces()[faceI].reverseFace();
        else face_ = mesh_.faces()[faceI];

    face_size_ = face_.size();

    pointStatus_.resize_nocopy(face_size_);

    // loop face
    forAll(face_, i)
    {        
        label ptIndex = pLabels.find(face_[i]);

        pointStatus_[i] = distanceList[ptIndex];
    }
    weight0_.resize_fill(face_size_, 0);
    weight1_.resize_fill(face_size_, 0);
}

void Foam::cutFacePLICNew::calcFaceNormal()
{
    faceNormal_ = (points_[face_[2]]-points_[face_[1]]) ^
                (points_[face_[0]]-points_[face_[1]]) ;
    
    // in case of a defect face
    if (mag(faceNormal_) == 0) faceNormal_ = face_.areaNormal(points_);

    if (mag(faceNormal_) == 0) faceNormal_ = Zero;
    else faceNormal_ = faceNormal_ / mag(faceNormal_);

}


void Foam::cutFacePLICNew::calcWeights(
    const scalar cutValue
)
{
    scalar phi_m, phi_m1;

    m1_found_ = -1;
    m2_found_ = -1;


    forAll(face_ , i)
    {
        label i1 = (i+1) % face_size_;
        label status= map_status(pointStatus_[i] - cutValue, pointStatus_[i1] - cutValue);
        
        switch (status)
            {
            case 1:
                break;
            case -1: 
                weight0_[i] = 1.0;
                break;
            case 0:
                phi_m = pointStatus_[i];
                phi_m1 = pointStatus_[i1];
                if (phi_m-cutValue < 0)
                {
                    weight0_[i] = phi_m / (phi_m - phi_m1);
                    weight1_[i] = 1.0 / (phi_m1 - phi_m);
                    if (m1_found_ == -1)
                    {
                        m1_found_ = i;                        
                        x0_k_const_ = points_[face_[i]]  + weight0_[i] * (points_[face_[i1]] - points_[face_[i]]  );
                        x0_k_lin_ = weight1_[i] * (points_[face_[i1]] - points_[face_[i]] );
                    }  
                }
                else
                {
                    weight0_[i] = phi_m1 / (phi_m1 - phi_m);

                    weight1_[i] = 1.0 / (phi_m - phi_m1);
                    if (m1_found_ == -1)
                    {
                        m1_found_ = i;
                        x0_k_const_ = points_[face_[i1]] + weight0_[i] * (points_[face_[i]]  - points_[face_[i1]] );
                        x0_k_lin_ = weight1_[i] * (points_[face_[i]]  - points_[face_[i1]] );
                    }
                }
                break;
            case -2:
                if (m2_found_ == -1)
                {
                    m2_found_ = i;
                    x0_k_lin_ = Zero;
                    if (pointStatus_[i]-cutValue > -eps_ && pointStatus_[i]-cutValue < eps_)
                        x0_k_const_ = points_[face_[i]] ;
                    else 
                        x0_k_const_ = points_[face_[i1]];
                }
                weight0_[i] = 1.0;
                weight1_[i] = 1.0/mag(pointStatus_[i]-pointStatus_[i1]); // or 0?
                break;
            case 2:
                if (m2_found_ == -1)
                {
                    m2_found_ = i;
                    x0_k_lin_ = Zero;
                    if (pointStatus_[i]-cutValue > -eps_ && pointStatus_[i]-cutValue < eps_)
                        x0_k_const_ = points_[face_[i]];
                    else
                        x0_k_const_ = points_[face_[i1]];
                }
                break;
            case 3:
                break;
            default:
                WarningInFunction << "no satus for edge found" <<endl;
            }
        
    }

}

Foam::label Foam::cutFacePLICNew::calcSubFace
(
    scalar cutValue
)
{    
    // scalar representing the edge normal in the 2d plane (after projection) which is cross product of edge and face normal
    scalar sc;
    // derivative of sc wrt cutValue
    scalar sc_d;

    scalar sum = 0, sum_d = 0, sum_dd = 0;

    // find unit vector to project on
    label pos = 0, l=1;
    scalar max = mag(faceNormal_[0]);
    while(l < 3){
        if (mag(faceNormal_[l]) > max){
            max = mag(faceNormal_[l]);
            pos = l;
        }
        l++;
    }
    label ind0 = (pos+1) % 3;
    label ind1 = (pos+2) % 3;

    // set x0_k for cases without intersecting edge
    if ((m1_found_ == -1)){   
        x0_k_lin_ = Zero;
        x0_k_const_ = points_[face_[0]];
        m1_found_ = 0;
        m2_found_ = face_size_ - 1;
    }

    // loop over edges to calculate area of the i-th face and its derivatives
    forAll (face_, i)
    {
        if (m1_found_ == i || m2_found_ == i)
            continue;  

        label i1 = (i+1) % face_size_;
        
        sc = (points_[face_[i]][ind0] - x0_k_const_[ind0] - cutValue * x0_k_lin_[ind0]) * (points_[face_[i1]][ind1] - points_[face_[i]][ind1]) +
            (points_[face_[i]][ind1] - x0_k_const_[ind1] - cutValue * x0_k_lin_[ind1]) * (points_[face_[i]][ind0] - points_[face_[i1]][ind0]);

        sc_d = -x0_k_lin_[ind0] * (points_[face_[i1]][ind1] - points_[face_[i]][ind1]) -
            x0_k_lin_[ind1] * (points_[face_[i]][ind0] - points_[face_[i1]][ind0]);

        sum += sc * (weight0_[i] + cutValue * weight1_[i]);

        sum_d += sc * weight1_[i] + sc_d * (weight0_[i] + cutValue * weight1_[i]);

        sum_dd += sc_d * weight1_[i];
    }
    
    subFaceArea_ = sum * 0.5 / faceNormal_[pos];

    subFaceAread_ = sum_d * 0.5 / faceNormal_[pos];
    subFaceAreadd_ = sum_dd / faceNormal_[pos];
        
    return 0;
}


// ************************************************************************* //

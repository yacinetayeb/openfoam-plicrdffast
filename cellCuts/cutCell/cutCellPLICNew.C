
#include "cutCellPLICNew.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutCellPLICNew::cutCellPLICNew(const fvMesh& mesh)
:
    cutCell(mesh),
    mesh_(mesh),
    points_(mesh.points()),
    tol_(1e-10),
    cellI_(-1),
    normal_(Zero),
    cutValue_(0),
    cutFace_(mesh_),
    faceCentre_(Zero),
    subCellVolume_(-10),
    VOF_(-10),
    k1_(-1)
{
    //clearStorage();
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void Foam::cutCellPLICNew::findFirstIntersecEdgeAndCalcX0
(
    const cell c,
    const scalar cutValue,
    const scalarList& distanceList
)
{

    // first and second point of the intersecting edge
    label pt1 =-1; label pt2 =-1;

    const labelList& pLabels = mesh_.cellPoints(cellI_);

    for (const label faceI : c)
    {
        // we go in this block after finding the first intersecting edge, to find k2 (the other face that has that same edge)
        if (k1_ != -1 )
            break;        

        // basically calcWeights, but reimplemented to break as soon as intersecting edge is found:
        scalar pointStatus_i1 = distanceList[pLabels.find(mesh_.faces()[faceI][0])] - cutValue;
        scalar eps=1e-12;
        forAll( mesh_.faces()[faceI], i )
        {
            label i1 = (i+1) % mesh_.faces()[faceI].size();

            scalar pointStatus_i = pointStatus_i1;
            pointStatus_i1 = distanceList[pLabels.find(mesh_.faces()[faceI][i1])] - cutValue;

            if (pointStatus_i < -eps && pointStatus_i1 > eps)
            {
                k1_ = faceI;
                pt1 = mesh_.faces()[faceI][i];
                pt2 = mesh_.faces()[faceI][i1];
                k1_m1_ = edge(pt1, pt2);

                x0_const_ = points_[pt1]  + (pointStatus_i + cutValue) / (pointStatus_i - pointStatus_i1) * (points_[pt2] - points_[pt1]);
                x0_lin_ = ( points_[pt2] - points_[pt1] ) / (pointStatus_i1 - pointStatus_i);
                break;
            }
            else if (pointStatus_i1 < -eps && pointStatus_i > eps)
            {
                k1_ = faceI;
                pt1 = mesh_.faces()[faceI][i];
                pt2 = mesh_.faces()[faceI][i1];
                k1_m1_ = edge(pt1, pt2);

                x0_const_ = points_[pt2] + (pointStatus_i1+cutValue) / (pointStatus_i1 - pointStatus_i) * (points_[pt1]  - points_[pt2] );
                x0_lin_ = (points_[pt1]  - points_[pt2] ) / (pointStatus_i - pointStatus_i1);
                break;
            }
            else if (pointStatus_i > -eps && pointStatus_i < eps)
            {
                k1_ = faceI;
                pt1 = mesh_.faces()[faceI][i];
                pt2 = mesh_.faces()[faceI][i1];
                k1_m1_ = edge(pt1, pt2);

                x0_const_ = points_[pt1] ;
                break;
            }
            else if (pointStatus_i1 > -eps && pointStatus_i1 < eps)
            {
                k1_ = faceI;
                pt1 = mesh_.faces()[faceI][i];
                pt2 = mesh_.faces()[faceI][i1];
                k1_m1_ = edge(pt1, pt2);

                x0_const_ = points_[pt2] ;
                break;
            }
        }
    
    }
}


// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cutCellPLICNew::calcSubCell
(
    const label celli,
    const scalar cutValue,
    const vector& normal,
    const point& base,
    const scalar volume, 
    const scalar tol,
    scalarList& distanceList
)
{
    //clearStorage();
    tol_ = tol;
    cellI_ = celli;
    cutValue_ = cutValue;
    normal_ = normal;

    const cell& c = mesh_.cells()[celli];
    //vector base = points_[mesh_.faces()[c[0]][0]];
    
    const labelList& pLabels = mesh_.cellPoints(celli);

    // only the case when called by mapAplhaField
    if (distanceList.empty()) {
        forAll(pLabels, pi)
        {
            scalar value = (points_[pLabels[pi]] - base) & normal;
            if (mag(value) < tol) value = 0;
            
            distanceList[pi] = value;
        }
    }

    faceCentre_ = base + normal_*cutValue_;

    k1_ = -1;
    
    findFirstIntersecEdgeAndCalcX0(c, cutValue, distanceList);

    // error handling for rare but possible cases
    if (k1_ == -1) {
        x0_lin_ = normal_;
        x0_const_ = base;
    }

    scalar sum_alpha = 0, sum_d_alpha = 0, sum_dd_alpha = 0, sum_ddd_alpha = 0;

    // loop over faces to calculate alpha and its derivatives
    for (const label faceI : c)
    {
        
        label found1 = mesh_.faces()[faceI].edgeDirection(k1_m1_);
        if (faceI == k1_  || found1)
            continue;  
        
        label owner=0;
            
        if (faceI < mesh_.owner().size())
            owner = mesh_.owner()[faceI];

        cutFace_.initCutFacePLICNew(cellI_, faceI, owner, distanceList, pLabels);

        cutFace_.calcWeights(cutValue_);

        cutFace_.calcFaceNormal();

        cutFace_.calcSubFace(cutValue_);

        scalar gamma_0 = (points_[cutFace_.getFace()[0]] - x0_const_) & cutFace_.faceNormal();
        scalar gamma_1 = -x0_lin_ & cutFace_.faceNormal();

        sum_alpha += (gamma_0 + cutValue_ * gamma_1) * cutFace_.subFaceArea();
        
        sum_d_alpha += (gamma_0 + cutValue_ * gamma_1) * cutFace_.subFaceAread() + gamma_1 * cutFace_.subFaceArea();

        sum_dd_alpha += (gamma_0 + cutValue_ * gamma_1) * cutFace_.subFaceAreadd() + 2 * gamma_1 * cutFace_.subFaceAread();

        sum_ddd_alpha += gamma_1 * cutFace_.subFaceAreadd();
    }

    subCellVolume_ = sum_alpha / 3.0;

    VOF_ = sum_alpha / 3.0 / volume;

    VOFd_ = sum_d_alpha / 3.0 / volume;

    VOFdd_ = sum_dd_alpha / 3.0 / volume;

    VOFddd_ = sum_ddd_alpha / volume;

    
    if (0 < mag(VOFddd_) && mag(VOFddd_) < tol_*tol_) VOFddd_ = 0;
    
    return 0;
}

void Foam::cutCellPLICNew::clearStorage()
{
    cellI_ = -1;
    cutValue_ = 0;
    x0_const_ = Zero;
    x0_lin_ = Zero;
    faceCentre_ = Zero;
    normal_ = Zero;
    subCellVolume_ = -10;
    VOF_ = -10;
}


// ************************************************************************* //

#ifndef cutFacePLICNew_H
#define cutFacePLICNew_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cutFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutFacePLICNew Declaration
\*---------------------------------------------------------------------------*/

class cutFacePLICNew
:
    public cutFace
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection
        //- with an isoface
        const fvMesh& mesh_;

        face face_;

        vector faceNormal_;

        scalar eps_;

        label face_size_;

        //- Storage for area vector of subface
        scalar subFaceArea_;

        //- Storage for area derivative vector of subface
        scalar subFaceAread_;

        //- Storage for area second derivative vector of subface
        scalar subFaceAreadd_;

        //- Storage for pointStatus_ cuts the cell at 0
        scalarList pointStatus_;

        //- Storage of the edge weight0
        scalarList weight0_;

        //- Storage of the edge weight1
        scalarList weight1_;

        //- Constant term of x0_k
        vector x0_k_const_;

        //- Linear term of x0_k
        vector x0_k_lin_;
        
        //- Storage for the first intersecting edge
        label m1_found_;

        //- Storage for m_1's neighbouring edge
        label m2_found_;

        const pointField& points_;

        // counter clockwise if true
        bool order_;


        //- A face status label taking one of the values:
        //
        //   -1: face is fully below the PLICsurface
        //    0: face is cut (has values larger and smaller than cutValue_)
        //   +1: face is fully above the PLICsurface
        label faceStatus_;
    
    // Private Member Functions

        label map_status(const scalar &status1, const scalar &status2);


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        explicit cutFacePLICNew(const fvMesh& mesh);


    // Member Functions
        //
        void initCutFacePLICNew(
            const label cellI,
            const label faceI,
            const label owner,
            const scalarList& distanceList,
            const labelList& pLabels
        );

        void calcFaceNormal();

        //
        void calcWeights(const scalar cutValue);

        //- Calculate cut points along edges of faceI
        label calcSubFace
        (
            scalar cutValue
        );

        const face& getFace() const noexcept
        {
            return face_;
        }

        const vector& faceNormal() const noexcept
        {
            return faceNormal_;
        }

        const label& m1_found() const noexcept
        {
            return m1_found_;
        }

        const label& m2_found() const noexcept
        {
            return m2_found_;
        }

        const vector& x0_k_const() noexcept
        {
            return x0_k_const_;
        }

        const vector& x0_k_lin() noexcept
        {
            return x0_k_lin_;
        }
        
        //- Returns area vector of cutted face
        const scalar& subFaceArea() const noexcept
        {
            return subFaceArea_;
        }
        
        //- Returns area derivative vector of cutted face
        const scalar& subFaceAread() const noexcept
        {
            return subFaceAread_;
        }
        
        //- Returns area second derivative vector of cutted face
        const scalar& subFaceAreadd() const noexcept
        {
            return subFaceAreadd_;
        }

        const scalarList& pointStatus() const noexcept
        {
            return pointStatus_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

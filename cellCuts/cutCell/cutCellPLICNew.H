#ifndef cutCellPLICNew_H
#define cutCellPLICNew_H

#include "cutCell.H"
#include "cutFacePLICNew.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutCellPLICNew Declaration
\*---------------------------------------------------------------------------*/

class cutCellPLICNew
:
    public cutCell
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;

        const pointField& points_;
        
        scalar tol_;

        //- Cell to cut
        label cellI_;

        //- Normal of the cutting plane, apperently looking in the fluid direction
        vector normal_;

        //- Cutvalue used to cut cell
        scalar cutValue_;

        //- An cutFacePLIC object to get access to its face cutting functionality
        cutFacePLICNew cutFace_;

        //- Face centre of the cutFace
        point faceCentre_;

        //- Face normal of the PLICface by convention pointing from high to low
        //  values (i.e. opposite of the gradient vector).
        vector faceArea_;

        //- Volume of fully submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VOF_;

        //- Volume of Fluid derivative for cellI
        scalar VOFd_;

        //- Volume of Fluid second derivative for cellI
        scalar VOFdd_;

        //- Volume of Fluid third derivative for cellI
        scalar VOFddd_;

        //- Constant term of x0
        vector x0_const_;

        //- Linear term of x0
        vector x0_lin_;

        // face with first intersecting edge
        label k1_;
        // first intersecting edge
        edge k1_m1_;
        

        //- A cell status label taking one of the values:
        //
        //   -1: cell is fully below the PLICsurface
        //    0: cell is cut
        //   +1: cell is fully above the PLICsurface
        label cellStatus_;

        void findFirstIntersecEdgeAndCalcX0
        (
            const cell c,
            const scalar cutValue,
            const scalarList& distanceList
        );


  public:

    // Constructors

        //- Construct from fvMesh
        explicit cutCellPLICNew(const fvMesh& mesh);


    // Member Functions

        //- Sets internal values and returns face status
        label calcSubCell
        (
            const label celli,
            const scalar cutValue,
            const vector& normal,
            const point& base,
            const scalar volume,
            const scalar tol,
            scalarList& distanceList
        );

        //- Returns subCellVolume
        scalar subCellVolume() const noexcept
        {
            return subCellVolume_;
        }

        //- Returns the centre of the cutting PLICface
        const point& faceCentre() const noexcept
        {
            return faceCentre_;
        }

        //- Returns the area normal vector of the cutting PLICface
        const vector& faceArea() const noexcept
        {
            return faceArea_;
        }

        //- sets the area normal vector of the cutting PLICface
        void setFaceArea(const vector& value) noexcept
        {
            faceArea_ = value;
        }

        //- Returns cellStatus per reference
        void setCellStatus(label value) noexcept
        {
            cellStatus_ = value;
        }

        //- Returns cellStatus
        label cellStatus() const noexcept
        {
            return cellStatus_;
        }

        //- Returns volume of fluid value
        scalar VolumeOfFluid() const noexcept
        {
            return VOF_;
        }
        
        //- Returns  volume of fluid derivative value
        scalar VolumeOfFluid_d() const noexcept
        {
            return VOFd_;
        }

        //- Returns volume of fluid second derivative value
        scalar VolumeOfFluid_dd() const noexcept
        {
            return VOFdd_;
        }

        //- Returns volume of fluid third derivative value
        scalar VolumeOfFluid_ddd() const noexcept
        {
            return VOFddd_;
        }

        //- Returns cutValue
        scalar cutValue() const noexcept
        {
            return cutValue_;
        }

        //- Resets internal values
        void clearStorage();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

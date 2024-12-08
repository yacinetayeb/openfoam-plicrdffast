/*---------------------------------------------------------------------------*\

Class
    Foam::surfaceIteratorPLICNew

Description
    Finds the cutValue that matches the volume fraction using Face-based PLIC
    positioning method from Bothe and Kromer

Author
    Yacine Tayeb, MagmaSoft

SourceFiles
    surfaceIteratorPLICNew.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceIteratorPLICNew_H
#define surfaceIteratorPLICNew_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cutCellPLICNew.H"
#include "cubicEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class surfaceIteratorPLICNew Declaration
\*---------------------------------------------------------------------------*/

class surfaceIteratorPLICNew
{
    // Private Data

        //- Mesh whose cells and faces to cut at their intersection with an
        //- isosurface.
        const fvMesh& mesh_;

        //- Cuts the cell and returns volume centre and normal
        cutCellPLICNew cutCell_;

        //- Tolerance for marking of surface cells:
        //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
        scalar surfCellTol_;

    // Private Member Functions

        //- Compute new step of the interface positioning algorithm 
        scalar compute_step(const scalar alpha1, const scalar tol);

        //- Cubic polynomial for the interface positioning algorithm
        scalar cubic_polynomial(const scalar s_i, const scalar s);


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        surfaceIteratorPLICNew(const fvMesh& mesh, const scalar tol);


    // Member Functions

        //- Determine if a cell is a surface cell
        bool isASurfaceCell(const scalar alpha1) const
        {
            return
            (
                surfCellTol_ < alpha1
             && alpha1 < 1 - surfCellTol_
            );
        }

        //- Finds matching cutValue for the given value fraction
        //  \return the cellStatus
        label vofCutCell
        (
            const label celli,
            const scalar alpha1,
            const scalar isoFaceTol,
            const label maxIter,
            vector normal
        );


        //- The volume of cutted volume
        scalar subCellVolume() const
        {
            return cutCell_.subCellVolume();
        }

        //- The centre of cutting isosurface
        const point& surfaceCentre() const
        {
            return cutCell_.faceCentre();
        }

        //- The area vector of cutting isosurface
        const vector& surfaceArea() const
        {
            return cutCell_.faceArea();
        }

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VolumeOfFluid() const
        {
            return cutCell_.VolumeOfFluid();
        }

        //- The cutValue
        scalar cutValue() const
        {
            return cutCell_.cutValue();
        }

        //- The cellStatus
        label cellStatus()
        {
            return cutCell_.cellStatus();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

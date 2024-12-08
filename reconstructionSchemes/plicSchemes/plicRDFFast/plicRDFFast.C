/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "plicRDFFast.H"
#include "interpolationCellPoint.H"
#include "fvc.H"
#include "leastSquareGrad.H"
#include "zoneDistribute.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(plicRDFFast, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,plicRDFFast, components);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reconstruction::plicRDFFast::gradSurf(const volScalarField& phi)
{
    addProfilingInFunction(geometricVoF);
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    exchangeFields.setUpCommforZone(interfaceCell_, false);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, mesh_.C())
    );
    Map<scalar> mapPhi
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);

    const labelListList& stencil = exchangeFields.getStencil();

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];

        cellCentre.clear();
        phiValues.clear();

        forAll(stencil[celli], j)
        {
            const label gblIdx = stencil[celli][j];
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields.getValue(phi, mapPhi, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];
        interfaceNormal_[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


void Foam::reconstruction::plicRDFFast::setInitNormals()
{
    addProfilingInFunction(geometricVoF);
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    interfaceLabels_.clear();

    forAll(alpha1_, celli)
    {
        if (sIterPLIC_.isASurfaceCell(alpha1_[celli]))
        {
            interfaceCell_[celli] = true; // is set to false earlier
            interfaceLabels_.append(celli);
        }
    }
    interfaceNormal_.setSize(interfaceLabels_.size());

    RDF_.markCellsNearSurf(interfaceCell_, 1);
    const boolList& nextToInterface_ = RDF_.nextToInterface();
    exchangeFields.updateStencil(nextToInterface_);

    gradSurf(alpha1_);
    
}

 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::plicRDFFast::plicRDFFast
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),

    interfaceNormal_(0.2*mesh_.nCells()),

    isoFaceTol_(modelDict().getOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().getOrDefault<scalar>("surfCellTol", 1e-8)),
    tol_(modelDict().getOrDefault<scalar>("tol", 1e-6)),
    relTol_(modelDict().getOrDefault<scalar>("relTol", 0.1)),
    iteration_(modelDict().getOrDefault<label>("iterations", 5)),
    RDF_(mesh_),
    sIterPLIC_(mesh_,surfCellTol_)
{
    setInitNormals();

    centre_ = dimensionedVector("centre", dimLength, Zero);
    normal_ = dimensionedVector("normal", dimArea, Zero);

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }
        
        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            20,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = Zero;
                centre_[celli] = Zero;
            }
        }
        else
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::plicRDFFast::reconstruct(bool forceUpdate)
{
    addProfilingInFunction(geometricVoF);
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    interfaceCell_ = false;

    // Sets interfaceCell_ and interfaceNormal
    setInitNormals();

    centre_ = dimensionedVector("centre", dimLength, Zero);
    normal_ = dimensionedVector("normal", dimArea, Zero);

    bitSet tooCoarse(mesh_.nCells(),false);

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(interfaceNormal_[i]) == 0 || tooCoarse.test(celli))
        {
            continue;
        }

        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            20,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = Zero;
                centre_[celli] = Zero;
            }
        }
        else
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
        }
    }

    normal_.correctBoundaryConditions();
    centre_.correctBoundaryConditions();
    
}


void Foam::reconstruction::plicRDFFast::mapAlphaField() const
{
    // without it we seem to get a race condition
    mesh_.C();

    cutCellPLICNew cutCell(mesh_);
    
    scalarList distanceList(0);

    Info << "mapAlphaField is called" << nl;

    forAll(normal_, celli)
    {
        if (mag(normal_[celli]) != 0)
        {
            vector n = normal_[celli]/mag(normal_[celli]);
            scalar cutValue = (centre_[celli] - mesh_.C()[celli]) & (n);
            cutCell.calcSubCell
            (
                celli,
                cutValue,
                n,
                mesh_.C()[celli],
                mesh_.V()[celli],
                isoFaceTol_,
                distanceList
            );
            alpha1_[celli] = cutCell.VolumeOfFluid();

        }
    }
    alpha1_.correctBoundaryConditions();
    alpha1_.oldTime () = alpha1_;
    alpha1_.oldTime().correctBoundaryConditions();
}


// ************************************************************************* //

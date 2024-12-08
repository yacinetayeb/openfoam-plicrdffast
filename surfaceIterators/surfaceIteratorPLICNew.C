
#include "surfaceIteratorPLICNew.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceIteratorPLICNew::surfaceIteratorPLICNew
(
    const fvMesh& mesh,
    const scalar tol
)
:
    mesh_(mesh),
    cutCell_(mesh_),
    surfCellTol_(tol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::surfaceIteratorPLICNew::vofCutCell
(
    const label celli,
    const scalar alpha1,
    const scalar isoFaceTol,
    const label maxIter,
    vector normal
)
{
    if (mag(normal) == 0)
    {
        WarningInFunction
            << "normal length is zero in cell: " << celli << nl
            << "try increasing nCorrectors" << endl;

        return sign(alpha1-0.5);
    }

    normal.normalise();

    normal = -normal;

    //Handling special case where method is handed an almost full or empty cell
    if (alpha1 < surfCellTol_)
    {
        cutCell_.setCellStatus(-1);
        return -1;
    }
    else if (1 - alpha1 < surfCellTol_)
    {
        cutCell_.setCellStatus(1);
        return 1;
    }

    cutCell_.setCellStatus(0);

    // Setting distance values from points to plane in order
    const labelList& pLabels = mesh_.cellPoints(celli);

    vector base  = mesh_.C()[celli];
    
    DynamicList<scalar> fvert;

    scalarList distanceList(pLabels.size());
    
    forAll(pLabels, pi)
    {
        scalar value = (mesh_.points()[pLabels[pi]] - base) & normal;
        if (mag(value) < isoFaceTol) value = 0;
        
        distanceList[pi] = value;
        // no duplicates
        fvert.push_uniq(value);
    }

    const labelList order(Foam::sortedOrder(fvert));

    scalar f1 = fvert[order.first()];
    scalar f2 = fvert[order.last()];

    //scalar volume = cutCell_.VolumeOfFluid();
    scalar volume = mesh_.V()[celli];

    // compute s_0
    
    scalar cutValue = f1 + (f2 - f1) * (0.5 - cos((acos(2 * alpha1 - 1) - 2 * M_PI) / 3));

    label n = 1;
    while (n < maxIter)
    {
        cutCell_.calcSubCell(celli, cutValue, normal, base, volume, surfCellTol_, distanceList);
        

        if (mag(cutCell_.VolumeOfFluid() - alpha1) < surfCellTol_)
            break;
        else
        {
            // this gives j the index s.t. s_j < cutValue <= s_j+1 bzw. si[j-1] <cutValue<= si[j]
            label j = 0;
            while (j < order.size() && fvert[order[j]] < cutValue)
                j++;

            if (j == 0 || j == fvert.size())
            {
                WarningInFunction << "j out of range! Make sure cutValue= " << cutValue << "isnt out of [s_,s+] = [" << f1 
                                << "," << f2 << "] for alpha = " << alpha1 << ", normal = " << normal 
                                << "and celli: " << celli << endl;
                break;
            }
            else
            {
                scalar alpha_minus = cubic_polynomial(fvert[order[j - 1]], cutValue);
                scalar alpha_plus = cubic_polynomial(fvert[order[j]], cutValue);

                if (mag(alpha_minus - alpha_plus) < SMALL) WarningInFunction << "alpha_ and alpha+ equal" << endl;
                
                // first assuming areas and derivs were just calculated with new cutValue, change later to lambda?

                if (alpha_minus <= alpha1 && alpha1 <= alpha_plus)
                {
                    cubicEqn S_n(cutCell_.VolumeOfFluid_ddd() / 6.0, cutCell_.VolumeOfFluid_dd() / 2.0, cutCell_.VolumeOfFluid_d(), cutCell_.VolumeOfFluid() - alpha1);
                    
                    Roots<3> roots =  S_n.roots();

                    for (scalar &root : roots)
                    {   
                        // check conversion from complex to double isnt causing problemo

                        root = root + cutValue; // + cutValue  because roots are for the polynom with variable X= cutValue* - cutValue
                        if (fvert[order[j-1]] - isoFaceTol <= root && root <= fvert[order[j]] + isoFaceTol)
                        {
                            cutValue = root;
                            
                            goto endofwhile;
                        }
                    }
                    WarningInFunction << "no root found!!!!!!!!!" << endl;

                } 
                else
                {
                    scalar delta_s = compute_step(alpha1, isoFaceTol);

                    if (f1 > cutValue + delta_s || cutValue + delta_s > f2) break;

                    cutValue += delta_s;
                    
                    n++;
                }
            }
        }

    }
    endofwhile:

    //return 0;
    
    // Check result
    label status = cutCell_.calcSubCell(celli, cutValue, normal, base, volume, surfCellTol_, distanceList);
    scalar VOF = cutCell_.VolumeOfFluid();

    scalar res = mag(VOF - alpha1);
    
    // area of enclosed plane segment is derivative of volume fraction wrt cutValue
    cutCell_.setFaceArea(-normal); //*cutCell_.VolumeOfFluid_d()    
    
    if (res > surfCellTol_)
    {

        WarningInFunction << "After " << n << " steps, VOF different from alpha1= " << alpha1 << " by " 
                        << res << " for normal: " << normal << endl
                                << "and celli: " << celli << endl;
        cutValue = f1 + (f2 - f1) * (0.5 - cos((acos(2 * alpha1 - 1) - 2 * M_PI) / 3));
        cutCell_.calcSubCell(celli, cutValue, normal, base, volume, surfCellTol_, distanceList);
        VOF = cutCell_.VolumeOfFluid();
        res = mag(VOF - alpha1);
        WarningInFunction << "taking first s0 guess: " << cutValue << " with VOF= " << VOF << " and error= " << res << endl; 

    }
    
    return 0;

    // If tolerance not met use the secant method  with f3 as a hopefully very
    // good initial guess to crank res the last piece down below tol

    /*scalar x2 = cutValue;
    scalar g2 = VOF - alpha1;
    scalar x1 = max(1e-3*(f2 - f1), 100*SMALL);
    x1 = max(x1, f1);
    x1 = min(x1, f2);
    cutCell_.calcSubCell(celli, x1,normal, base, volume, tol);
    scalar g1 = cutCell_.VolumeOfFluid() - alpha1;

    n = 0;
    scalar g0(0), x0(0);
    while (res > tol && n < maxIter && mag(g1 - g2) > tol)
    {
        x0 = (x2*g1 - x1*g2)/(g1 - g2);
        status = cutCell_.calcSubCell(celli, x0, normal, base, volume, tol);
        g0 = cutCell_.VolumeOfFluid() - alpha1;
        res = mag(g0);
        x2 = x1; g2 = g1;
        x1 = x0; g1 = g0;
        n++;
    }

    return status;*/
}


Foam::scalar Foam::surfaceIteratorPLICNew::compute_step(
    const scalar alpha1,
    const scalar tol)
{

    scalar discriminant = cutCell_.VolumeOfFluid_d() * cutCell_.VolumeOfFluid_d() - 2 * (cutCell_.VolumeOfFluid() - alpha1) * cutCell_.VolumeOfFluid_dd();
    scalar delta_s;
    if (discriminant >= 0 && mag(cutCell_.VolumeOfFluid_dd()) >= tol)
    {
        delta_s = (sqrt(discriminant) - cutCell_.VolumeOfFluid_d()) / cutCell_.VolumeOfFluid_dd();
        
    }
    else
    {
        if (mag(cutCell_.VolumeOfFluid_d()) > SMALL)
            delta_s = (alpha1 - cutCell_.VolumeOfFluid()) / cutCell_.VolumeOfFluid_d();
        else 
        {
            WarningInFunction << "VOFd is zero while trying to calculate new step" << endl;

            // a small delta_s to help get away from the cutValue where VOF is "constant"
            delta_s = 2*tol;
        }
    }
    return delta_s;
}

Foam::scalar Foam::surfaceIteratorPLICNew::cubic_polynomial(const scalar s_i, const scalar cutValue)
{
    scalar x = s_i - cutValue;
    scalar x2 = x * x;
    scalar x3 = x2 * x;
    return cutCell_.VolumeOfFluid_ddd() / 6.0 * x3 + cutCell_.VolumeOfFluid_dd() / 2.0 * x2 + cutCell_.VolumeOfFluid_d() * x + cutCell_.VolumeOfFluid();
}

// ************************************************************************* //

/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef dkmshell_h
#define dkmshell_h

#include "sm/Elements/structuralelement.h"
#include "floatmatrixf.h"

#define _IFT_DKMTShell_Name "dkmtshell"
#define _IFT_DKMQShell_Name "dkmqshell"
#define _IFT_DKMShell_kgw "kgw"
#define _IFT_DKMShell_kirchhoff "kirchhoff"

namespace oofem {
class ParamKey;

/**
 * Flat shell element with 3 (DKMTShell) or 4 (DKMQShell) nodes and 6 DOFs per node.
 * The membrane is the Allman / Ibrahimbegovic-Taylor-Wilson (ITW) element with a true drilling
 * rotation; the plate is the discrete Kirchhoff-Mindlin element of Katili (DKMT, DKMQ), free of
 * shear locking. The initial stress matrix is provided, so the element can be used in linear
 * stability and second order analyses.
 *
 * Matrices and stress resultants are computed by the plain C kernel in dkmshellkernel.c.
 * Only linear elastic isotropic material (IsoLE) and SimpleCS with thickness are supported.
 *
 * Optional input fields:
 * - kgw 1: gradient of w in the initial stress matrix from the nodal w only (default: from the rotations);
 * - kirchhoff 1: transverse shear deformation neglected (the plate becomes DKT / DKQ), for comparisons.
 *
 * Local system: e3 is the normal oriented by the node numbering (counterclockwise seen from +e3);
 * e1 is the side 1->2 projected on the mean plane. Warped quadrilaterals are projected on the mean
 * plane; the largest distance of a node from that plane is printed.
 *
 * Output (printOutputAt), per element: the local axes and, for every vertex, the stress resultants
 * per unit width in the element local system, evaluated directly from the element field:
 *     SHELL axes e1x e1y e1z  e2x e2y e2z  e3x e3y e3z
 *     SHELL vertex <k> node <label> Nx Ny Nxy Mx My Mxy Tx Ty
 *     SHELL warp <distance>                  (non-planar quadrilaterals only)
 */
class DKMShell : public StructuralElement
{
protected:
    int nn;
    /// Option of the initial stress matrix (DKMSHELL_KG_W_ROTATIONS or DKMSHELL_KG_W_LINEAR)
    int kgOption = 0;
    /// 1: transverse shear deformation neglected (DKT/DKQ)
    int kirchhoff = 0;
    /// Rows: e1, e2, e3 (global -> local)
    FloatMatrixF< 3, 3 >R;
    /// Local node coordinates in the mean plane
    double xl [ 4 ], yl [ 4 ];
    /// Largest distance of the nodes from the mean plane (quadrilaterals)
    double warp = 0.0;
    bool geometryReady = false;

    static ParamKey IPK_DKMShell_kgw;
    static ParamKey IPK_DKMShell_kirchhoff;

    void computeLocalGeometry();
    void giveMaterialData(double &E, double &nu, double &t, double &rho);

public:
    DKMShell(int n, Domain *d, int nnode);
    virtual ~DKMShell() { }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;
    int computeNumberOfDofs() override { return 6 * nn; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode) override;

    /// Stress resultants at the vertices in the local system: S[8*i + c] (see dkmshellkernel.h)
    void computeVertexStresses(TimeStep *tStep, double *S);
    /// Local axes e1, e2, e3 (9 values, by rows)
    void giveLocalAxes(double *e);
    void printOutputAt(FILE *file, TimeStep *tStep) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    double computeVolumeAround(GaussPoint *gp) override;
    double giveArea();
    integrationDomain giveIntegrationDomain() const override { return nn == 3 ? _Triangle : _Square; }
    Element_Geometry_Type giveGeometryType() const override { return nn == 3 ? EGT_triangle_1 : EGT_quad_1; }
    MaterialMode giveMaterialMode() override { return _PlaneStress; }
    int testElementExtension(ElementExtension ext) override { return 0; }

protected:
    void computeGaussPoints() override;
    // Not used: stiffness, mass and stress resultants come from the kernel
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS) override { answer.clear(); }
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { answer.clear(); }
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override { answer.clear(); }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override { answer.clear(); }
};

/// Three-node flat shell: ITW membrane + DKMT plate
class DKMTShell : public DKMShell
{
public:
    DKMTShell(int n, Domain *d) : DKMShell(n, d, 3) { }
    const char *giveInputRecordName() const override { return _IFT_DKMTShell_Name; }
    const char *giveClassName() const override { return "DKMTShell"; }
};

/// Four-node flat shell: ITW membrane + DKMQ plate
class DKMQShell : public DKMShell
{
public:
    DKMQShell(int n, Domain *d) : DKMShell(n, d, 4) { }
    const char *giveInputRecordName() const override { return _IFT_DKMQShell_Name; }
    const char *giveClassName() const override { return "DKMQShell"; }
};
} // end namespace oofem
#endif // dkmshell_h

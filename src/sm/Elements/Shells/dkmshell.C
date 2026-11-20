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

#include "sm/Elements/Shells/dkmshell.h"
#include "sm/Elements/Shells/dkmshellkernel.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "node.h"
#include "load.h"
#include "classfactory.h"
#include "crosssection.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "dofmanager.h"
#include "error.h"
#include "parametermanager.h"
#include "paramkey.h"
#include "domain.h"
#include <cmath>

namespace oofem {
REGISTER_Element(DKMTShell);
REGISTER_Element(DKMQShell);
ParamKey DKMShell::IPK_DKMShell_kgw(_IFT_DKMShell_kgw);
ParamKey DKMShell::IPK_DKMShell_kirchhoff(_IFT_DKMShell_kirchhoff);

DKMShell :: DKMShell(int n, Domain *d, int nnode) : StructuralElement(n, d), nn(nnode)
{
    numberOfDofMans = nnode;
    numberOfGaussPoints = 1;
}

void
DKMShell :: initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority)
{
    StructuralElement :: initializeFrom(ir, priority);
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    PM_UPDATE_PARAMETER(kgOption, ppm, ir, this->number, IPK_DKMShell_kgw, priority);
    PM_UPDATE_PARAMETER(kirchhoff, ppm, ir, this->number, IPK_DKMShell_kirchhoff, priority);
}

void
DKMShell :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u, D_v, D_w, R_u, R_v, R_w };
}

void
DKMShell :: computeGaussPoints()
{
    // A single point at the centroid: only needed by the general OOFEM machinery (material status, output)
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std :: make_unique< GaussIntegrationRule >(1, this, 1, 3);
        if ( nn == 3 ) {
            integrationRulesArray [ 0 ]->SetUpPointsOnTriangle(1, _PlaneStress);
        } else {
            integrationRulesArray [ 0 ]->SetUpPointsOnSquare(1, _PlaneStress);
        }
    }
}

void
DKMShell :: computeLocalGeometry()
{
    if ( geometryReady ) {
        return;
    }
    FloatArrayF< 3 >X [ 4 ], c;
    for ( int i = 0; i < nn; i++ ) {
        X [ i ] = FloatArrayF< 3 >(this->giveNode(i + 1)->giveCoordinates() );
        c += X [ i ] * ( 1.0 / nn );
    }
    FloatArrayF< 3 >e1, e2, e3;
    if ( nn == 3 ) {
        e1 = X [ 1 ] - X [ 0 ];
        e3 = cross(e1, X [ 2 ] - X [ 0 ]);
    } else {
        e3 = cross(X [ 2 ] - X [ 0 ], X [ 3 ] - X [ 1 ]);
        e1 = X [ 1 ] - X [ 0 ];
    }
    double n3 = norm(e3);
    if ( n3 <= 0.0 ) {
        OOFEM_ERROR("element %d: degenerate geometry", this->giveLabel() );
    }
    e3 *= 1.0 / n3;
    e1 = e1 - e3 * dot(e1, e3);
    e1 *= 1.0 / norm(e1);
    e2 = cross(e3, e1);
    for ( int j = 0; j < 3; j++ ) {
        R(0, j) = e1 [ j ];
        R(1, j) = e2 [ j ];
        R(2, j) = e3 [ j ];
    }
    warp = 0.0;
    for ( int i = 0; i < nn; i++ ) {
        FloatArrayF< 3 >d = X [ i ] - c;
        xl [ i ] = dot(d, e1);
        yl [ i ] = dot(d, e2);
        warp = std :: max(warp, fabs(dot(d, e3) ) );
    }
    geometryReady = true;
}

void
DKMShell :: giveLocalAxes(double *e)
{
    this->computeLocalGeometry();
    for ( int r = 0; r < 3; r++ ) {
        for ( int c = 0; c < 3; c++ ) {
            e [ 3 * r + c ] = R(r, c);
        }
    }
}

bool
DKMShell :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    this->computeLocalGeometry();
    answer.resize(6 * nn, 6 * nn);
    answer.zero();
    for ( int i = 0; i < nn; i++ ) {
        for ( int b = 0; b < 2; b++ ) {
            for ( int r = 0; r < 3; r++ ) {
                for ( int s = 0; s < 3; s++ ) {
                    answer.at(6 * i + 3 * b + r + 1, 6 * i + 3 * b + s + 1) = R(r, s);
                }
            }
        }
    }
    return true;
}

int
DKMShell :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
{
    // body loads (3 global force components)
    this->computeLocalGeometry();
    answer.resize(3, 3);
    for ( int r = 0; r < 3; r++ ) {
        for ( int s = 0; s < 3; s++ ) {
            answer.at(r + 1, s + 1) = R(r, s);
        }
    }
    return 1;
}

void
DKMShell :: giveMaterialData(double &E, double &nu, double &t, double &rho)
{
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    auto *mat = dynamic_cast< IsotropicLinearElasticMaterial * >( this->giveStructuralCrossSection()->giveMaterial(gp) );
    if ( !mat ) {
        OOFEM_ERROR("element %d: isotropic linear elastic material (IsoLE) required", this->giveLabel() );
    }
    E = mat->giveYoungsModulus();
    nu = mat->givePoissonsRatio();
    t = this->giveCrossSection()->give(CS_Thickness, gp);
    rho = this->giveCrossSection()->give('d', gp);
}

static void setKernelMaterial(DKMShellMaterial &m, double E, double nu, double t, int kirchhoff = 0)
{
    m.E = E;
    m.nu = nu;
    m.t = t;
    m.k = kirchhoff ? 1.0e8 : 5.0 / 6.0;    // huge k: phi -> 0, the plate becomes DKT/DKQ
    m.gdrill = 0.0;
}

void
DKMShell :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    this->computeLocalGeometry();
    double E, nu, t, rho;
    this->giveMaterialData(E, nu, t, rho);
    DKMShellMaterial m;
    setKernelMaterial(m, E, nu, t, kirchhoff);
    int nd = 6 * nn;
    double K [ 24 * 24 ];
    if ( dkmshell_stiffness(nn, xl, yl, & m, K) ) {
        OOFEM_ERROR("element %d: invalid geometry (clockwise nodes or distorted element)", this->giveLabel() );
    }
    answer.resize(nd, nd);
    for ( int i = 0; i < nd; i++ ) {
        for ( int j = 0; j < nd; j++ ) {
            answer.at(i + 1, j + 1) = K [ i * nd + j ];
        }
    }
}

void
DKMShell :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    this->computeLocalGeometry();
    double E, nu, t, rho;
    this->giveMaterialData(E, nu, t, rho);
    DKMShellMaterial m;
    setKernelMaterial(m, E, nu, t, kirchhoff);
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);      // already in the local system
    int nd = 6 * nn;
    double ul [ 24 ], KG [ 24 * 24 ];
    for ( int i = 0; i < nd; i++ ) {
        ul [ i ] = u.at(i + 1);
    }
    if ( dkmshell_initial_stress(nn, xl, yl, & m, ul, kgOption, KG) ) {
        OOFEM_ERROR("element %d: invalid geometry", this->giveLabel() );
    }
    answer.resize(nd, nd);
    for ( int i = 0; i < nd; i++ ) {
        for ( int j = 0; j < nd; j++ ) {
            answer.at(i + 1, j + 1) = KG [ i * nd + j ];
        }
    }
}

double
DKMShell :: giveArea()
{
    this->computeLocalGeometry();
    return dkmshell_area(nn, xl, yl);
}

double
DKMShell :: computeVolumeAround(GaussPoint *gp)
{
    return this->giveArea() / this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}

void
DKMShell :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // mass lumped on the translations: rho t A / n per node (no rotational mass)
    double E, nu, t, rho;
    this->giveMaterialData(E, nu, t, rho);
    double mn = rho * t * this->giveArea() / nn;
    answer.resize(6 * nn, 6 * nn);
    answer.zero();
    for ( int i = 0; i < nn; i++ ) {
        for ( int c = 1; c <= 3; c++ ) {
            answer.at(6 * i + c, 6 * i + c) = mn;
        }
    }
}

void
DKMShell :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    this->computeLumpedMassMatrix(answer, tStep);
}

void
DKMShell :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    // dead weight and similar: force per unit mass, global -> local, split equally among the nodes
    FloatArray force;
    FloatMatrix T;
    forLoad->computeComponentArrayAt(force, tStep, mode);
    answer.clear();
    if ( force.giveSize() < 3 ) {
        return;
    }
    FloatArray f3 = FloatArray::fromIniList({ force.at(1), force.at(2), force.at(3) });
    this->computeLoadGToLRotationMtrx(T);
    f3.rotatedWith(T, 'n');
    double E, nu, t, rho;
    this->giveMaterialData(E, nu, t, rho);
    double mn = rho * t * this->giveArea() / nn;
    answer.resize(6 * nn);
    answer.zero();
    for ( int i = 0; i < nn; i++ ) {
        for ( int c = 1; c <= 3; c++ ) {
            answer.at(6 * i + c) = mn * f3.at(c);
        }
    }
}

void
DKMShell :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // linear elastic: f = K u in the local system
    FloatMatrix K;
    FloatArray u;
    this->computeStiffnessMatrix(K, TangentStiffness, tStep);
    this->computeVectorOf(VM_Total, tStep, u);
    answer.beProductOf(K, u);
}

void
DKMShell :: computeVertexStresses(TimeStep *tStep, double *S)
{
    this->computeLocalGeometry();
    double E, nu, t, rho;
    this->giveMaterialData(E, nu, t, rho);
    DKMShellMaterial m;
    setKernelMaterial(m, E, nu, t, kirchhoff);
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    double ul [ 24 ];
    for ( int i = 0; i < 6 * nn; i++ ) {
        ul [ i ] = u.at(i + 1);
    }
    dkmshell_vertex_stresses(nn, xl, yl, & m, ul, S);
}

void
DKMShell :: printOutputAt(FILE *file, TimeStep *tStep)
{
    fprintf(file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
    if ( !this->isActivated(tStep) ) {
        fprintf(file, "is not active in current time step\n");
        return;
    }
    this->computeLocalGeometry();
    fprintf(file, "  SHELL axes %.10e %.10e %.10e  %.10e %.10e %.10e  %.10e %.10e %.10e\n",
            R(0, 0), R(0, 1), R(0, 2), R(1, 0), R(1, 1), R(1, 2), R(2, 0), R(2, 1), R(2, 2) );
    double S [ 4 * DKMSHELL_NS ];
    this->computeVertexStresses(tStep, S);
    for ( int i = 0; i < nn; i++ ) {
        fprintf(file, "  SHELL vertex %d node %d", i + 1, this->giveNode(i + 1)->giveLabel() );
        for ( int c = 0; c < DKMSHELL_NS; c++ ) {
            fprintf(file, " %.10e", S [ DKMSHELL_NS * i + c ]);
        }
        fprintf(file, "\n");
    }
    if ( warp > 0.0 ) {
        fprintf(file, "  SHELL warp %.6e\n", warp);
    }
}

int
DKMShell :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    // values at the centroid, in the local system, in the OOFEM shell tensor format
    if ( type == IST_ShellForceTensor || type == IST_ShellMomentTensor ) {
        this->computeLocalGeometry();
        double E, nu, t, rho;
        this->giveMaterialData(E, nu, t, rho);
        DKMShellMaterial m;
        setKernelMaterial(m, E, nu, t, kirchhoff);
        FloatArray u;
        this->computeVectorOf(VM_Total, tStep, u);
        double ul [ 24 ], S [ DKMSHELL_NS ];
        for ( int i = 0; i < 6 * nn; i++ ) {
            ul [ i ] = u.at(i + 1);
        }
        double c = nn == 3 ? 1.0 / 3.0 : 0.0;
        dkmshell_point_stresses(nn, xl, yl, & m, ul, c, c, S);
        if ( type == IST_ShellForceTensor ) {
            answer = FloatArray::fromIniList({ S [ 0 ], S [ 1 ], 0.0, S [ 7 ], S [ 6 ], S [ 2 ] });
        } else {
            answer = FloatArray::fromIniList({ S [ 3 ], S [ 4 ], 0.0, 0.0, 0.0, S [ 5 ] });
        }
        return 1;
    }
    return StructuralElement :: giveIPValue(answer, gp, type, tStep);
}
} // end namespace oofem

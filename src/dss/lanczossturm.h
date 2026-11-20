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

#ifndef lanczossturm_h
#define lanczossturm_h

#include "sparsegeneigenvalsystemnm.h"
#include "convergedreason.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "floatmatrix.h"

#define _IFT_LanczosSturm_Name "lanczossturm"

namespace oofem {
class Domain;
class EngngModel;

/**
 * Generalized eigenvalue solver K y = lambda M y.
 * Shift-invert Lanczos (operator (K - sigma M)^-1 M) with full reorthogonalization in the M inner
 * product, Sturm sequence check (number of negative pivots of K - sigma_c M) and restart with
 * deflation until the number of eigenvalues found below sigma_c equals the Sturm count.
 * Eigenvectors are purified when M is singular (Nour-Omid, Parlett, Ericsson, Jensen 1987).
 * Supported matrices: DSSMatrix (smtype 8) or Skyline (smtype 0).
 *
 * Selected with stype 3 (GES_LanczosSturm). Optional environment variables:
 * - OOFEM_LANCZOS_SHIFT: initial shift (default 0);
 * - OOFEM_LANCZOS_MAXSTEPS: maximum number of Lanczos steps per round (default max(6 nroot + 60, 300)).
 */
class OOFEM_EXPORT LanczosSturm : public SparseGeneralEigenValueSystemNM
{
public:
    LanczosSturm(Domain *d, EngngModel *m);
    virtual ~LanczosSturm() { }

    ConvergedReason solve(SparseMtrx &A, SparseMtrx &B, FloatArray &x, FloatMatrix &v, double rtol, int nroot) override;
    const char *giveClassName() const override { return "LanczosSturm"; }
};
} // end namespace oofem
#endif // lanczossturm_h

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

#ifndef responsespectrum_h
#define responsespectrum_h

#include "sm/EngineeringModels/structengngmodel.h"
#include "geneigvalsolvertype.h"
#include "sparselinsystemnm.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrxtype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include <vector>

///@name Input fields for ResponseSpectrum
//@{
#define _IFT_ResponseSpectrum_Name "responsespectrum"
#define _IFT_ResponseSpectrum_nroot "nroot"
#define _IFT_ResponseSpectrum_rtolv "rtolv"
#define _IFT_ResponseSpectrum_periods "spectrumperiods"
#define _IFT_ResponseSpectrum_values "spectrumvalues"
#define _IFT_ResponseSpectrum_scale "spectrumscale"
#define _IFT_ResponseSpectrum_directions "directions"
#define _IFT_ResponseSpectrum_damping "damping"
#define _IFT_ResponseSpectrum_combination "modalcombination"
#define _IFT_ResponseSpectrum_minMassRatio "minmassratio"
#define _IFT_ResponseSpectrum_maxNroot "maxnroot"
#define _IFT_ResponseSpectrum_pdelta "pdelta"
//@}

namespace oofem {
/**
 * Linear dynamic analysis with a response spectrum (modal response spectrum analysis,
 * e.g. Eurocode 8 4.3.3.3, Italian NTC 2018 7.3.3.1).
 *
 * Input record:
 *     ResponseSpectrum nroot <n> [rtolv <tol>] [smtype <s>] [lstype <l>]
 *         spectrumperiods <n> T1 T2 ...   spectrumvalues <n> Sa1 Sa2 ...   [spectrumscale <f>]
 *         [directions <n> d1 d2 ...]      1 = X, 2 = Y, 3 = Z (default 2 1 2)
 *         [damping <xi>]                  viscous damping ratio for CQC, default 0.05
 *         [modalcombination <0|1>]        0 = CQC (default), 1 = SRSS
 *         [minmassratio <f>]              default 0.85: the number of modes is doubled until the cumulative
 *                                         effective mass reaches f in every requested direction
 *         [maxnroot <n>]                  upper limit of the number of modes, default 400
 *         [pdelta]                        second order: modes of K + KG, where KG comes from a static solution
 *                                         of all loads in the input (the gravity loads of the seismic
 *                                         combination, evaluated at t = 1)
 *
 * Masses come from the elements (density) and from lumpedmass elements. The spectrum Sa(T) is an
 * acceleration in units consistent with the model, interpolated linearly and constant outside the given
 * range. Eigenpairs are computed by the Lanczos solver with Sturm sequence check (GES_LanczosSturm).
 *
 * Output (.out file): table of modes with participation factors and effective masses, then for every
 * direction the combined values, with the sign of the largest modal contribution:
 *     RSA direction <X|Y|Z>
 *     RSA base_shear <V>
 *     RSA node <label> u v w rx ry rz
 *     RSA forces <label> <n> f1 ... fn          element end forces in the local system (K u)
 *     RSA shell <label> axes e1 e2 e3            DKMTShell / DKMQShell, then one line per vertex:
 *     RSA vertex <k> node <label> Nx Ny Nxy Mx My Mxy Tx Ty
 */
class ResponseSpectrum : public StructuralEngngModel
{
protected:
    int nroot = 12, maxNroot = 400, combination = 0;
    double rtolv = 1e-8, damping = 0.05, minMassRatio = 0.85, scale = 1.0;
    bool pdelta = false;
    FloatArray spectrumT, spectrumSa;
    IntArray directions;
    SparseMtrxType sparseMtrxType = SMT_Skyline;
    LinSystSolverType solverType = ST_Direct;

    /// Eigenvalues
    FloatArray omega2;
    /// Mass-normalized eigenvectors (neq x nmodes)
    FloatMatrix eigVec;
    /// Participation factors (nmodes x 3)
    FloatMatrix gamma;
    /// Total mass per direction (3)
    FloatArray totalMass;
    /// Static state of the gravity loads (pdelta)
    FloatArray uStatic;
    /// Unknowns served by giveUnknownComponent: 0 = static state, i > 0 = mode i
    int activeMode = 0;

    double giveSpectralAcceleration(double T) const;
    double giveCorrelationCoefficient(int i, int j) const;
    void assembleStiffness(SparseMtrx &K, TimeStep *tStep);
    void computeModes(TimeStep *tStep, int n);
    double combineModalValues(const std :: vector< double > &q, int &sign) const;

public:
    ResponseSpectrum(int i, EngngModel *master = nullptr);
    virtual ~ResponseSpectrum() { }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void solveYourself() override;
    void doStepOutput(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override { }
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    TimeStep *giveNextStep() override;
    bool newDofHandling() override { return true; }

    const char *giveClassName() const override { return "ResponseSpectrum"; }
    const char *giveInputRecordName() const { return _IFT_ResponseSpectrum_Name; }
};
} // end namespace oofem
#endif // responsespectrum_h

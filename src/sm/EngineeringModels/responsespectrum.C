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

#include "sm/EngineeringModels/responsespectrum.h"
#include "sm/Elements/Shells/dkmshell.h"
#include "sm/Elements/Shells/dkmshellkernel.h"
#include "sm/Elements/structuralelement.h"
#include "node.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "domain.h"
#include "classfactory.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "assemblercallback.h"
#include "oofemcfg.h"
#include "mathfem.h"
#include <cmath>
#include <cstdlib>
#include <ctime>

namespace oofem {
REGISTER_EngngModel(ResponseSpectrum);

ResponseSpectrum :: ResponseSpectrum(int i, EngngModel *master) : StructuralEngngModel(i, master)
{
    numberOfSteps = 1;
    ndomains = 1;
}

void
ResponseSpectrum :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    IR_GIVE_FIELD(ir, nroot, _IFT_ResponseSpectrum_nroot);
    IR_GIVE_OPTIONAL_FIELD(ir, rtolv, _IFT_ResponseSpectrum_rtolv);
    IR_GIVE_FIELD(ir, spectrumT, _IFT_ResponseSpectrum_periods);
    IR_GIVE_FIELD(ir, spectrumSa, _IFT_ResponseSpectrum_values);
    if ( spectrumT.giveSize() < 1 || spectrumT.giveSize() != spectrumSa.giveSize() ) {
        OOFEM_ERROR("spectrumperiods and spectrumvalues must have the same size");
    }
    IR_GIVE_OPTIONAL_FIELD(ir, scale, _IFT_ResponseSpectrum_scale);
    directions = { 1, 2 };
    IR_GIVE_OPTIONAL_FIELD(ir, directions, _IFT_ResponseSpectrum_directions);
    IR_GIVE_OPTIONAL_FIELD(ir, damping, _IFT_ResponseSpectrum_damping);
    IR_GIVE_OPTIONAL_FIELD(ir, combination, _IFT_ResponseSpectrum_combination);
    IR_GIVE_OPTIONAL_FIELD(ir, minMassRatio, _IFT_ResponseSpectrum_minMassRatio);
    IR_GIVE_OPTIONAL_FIELD(ir, maxNroot, _IFT_ResponseSpectrum_maxNroot);
    pdelta = ir->hasField(_IFT_ResponseSpectrum_pdelta);

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;
    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    if ( ( outputStream = fopen(this->dataOutputFileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("Can't open output file %s", this->dataOutputFileName.c_str() );
    }
    fprintf(outputStream, "%s", PRG_HEADER);
    fprintf(outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
    fprintf(outputStream, "%s\n", simulationDescription.c_str() );
}

TimeStep *
ResponseSpectrum :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;
    if ( currentStep ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }
    previousStep = std :: move(currentStep);
    // time 1: loads (gravity loads for P-Delta) are evaluated at t = 1
    currentStep = std :: make_unique< TimeStep >(istep, this, 1, 1.0, 1.0, counter);
    return currentStep.get();
}

double
ResponseSpectrum :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    int eq = dof->__giveEquationNumber();
    if ( eq <= 0 ) {
        return 0.0;                             // constrained DOF (zero prescribed value)
    }
    if ( activeMode > 0 ) {
        return eigVec.at(eq, activeMode);       // shape of the active mode
    }
    return uStatic.isNotEmpty() ? uStatic.at(eq) : 0.0;
}

double
ResponseSpectrum :: giveSpectralAcceleration(double T) const
{
    int n = spectrumT.giveSize();
    if ( T <= spectrumT.at(1) ) {
        return scale * spectrumSa.at(1);
    }
    for ( int i = 2; i <= n; i++ ) {
        if ( T <= spectrumT.at(i) ) {
            double a = ( T - spectrumT.at(i - 1) ) / ( spectrumT.at(i) - spectrumT.at(i - 1) );
            return scale * ( spectrumSa.at(i - 1) + a * ( spectrumSa.at(i) - spectrumSa.at(i - 1) ) );
        }
    }
    return scale * spectrumSa.at(n);
}

double
ResponseSpectrum :: giveCorrelationCoefficient(int i, int j) const
{
    // CQC correlation coefficient (Der Kiureghian), same damping for all modes
    if ( combination == 1 ) {
        return i == j ? 1.0 : 0.0;
    }
    double r = sqrt(omega2.at(j) / omega2.at(i) );
    double xi = damping;
    double num = 8.0 * xi * xi * ( 1.0 + r ) * pow(r, 1.5);
    double den = ( 1.0 - r * r ) * ( 1.0 - r * r ) + 4.0 * xi * xi * r * ( 1.0 + r ) * ( 1.0 + r );
    return num / den;
}

double
ResponseSpectrum :: combineModalValues(const std :: vector< double > &q, int &sign) const
{
    int n = ( int ) q.size();
    double s = 0.0, qmax = 0.0;
    sign = 1;
    for ( int i = 0; i < n; i++ ) {
        if ( fabs(q [ i ]) > qmax ) {
            qmax = fabs(q [ i ]);
            sign = q [ i ] >= 0 ? 1 : -1;
        }
        s += q [ i ] * q [ i ];
        if ( combination == 0 ) {
            for ( int j = i + 1; j < n; j++ ) {
                s += 2.0 * this->giveCorrelationCoefficient(i + 1, j + 1) * q [ i ] * q [ j ];
            }
        }
    }
    return sqrt(std :: max(s, 0.0) );
}

void
ResponseSpectrum :: assembleStiffness(SparseMtrx &K, TimeStep *tStep)
{
    K.buildInternalStructure(this, 1, EModelDefaultEquationNumbering() );
    this->assemble(K, tStep, TangentAssembler(TangentStiffness), EModelDefaultEquationNumbering(), this->giveDomain(1) );
    if ( pdelta ) {
        // second order: initial stress matrix of the gravity load state, added to K
        activeMode = 0;
        this->assemble(K, tStep, InitialStressMatrixAssembler(), EModelDefaultEquationNumbering(), this->giveDomain(1) );
    }
}

void
ResponseSpectrum :: computeModes(TimeStep *tStep, int n)
{
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering() );
    auto K = classFactory.createSparseMtrx(sparseMtrxType);
    auto M = classFactory.createSparseMtrx(sparseMtrxType);
    this->assembleStiffness(* K, tStep);
    M->buildInternalStructure(this, 1, EModelDefaultEquationNumbering() );
    this->assemble(* M, tStep, MassMatrixAssembler(), EModelDefaultEquationNumbering(), d);

    // Influence vectors (unit rigid body motion in X, Y, Z) and their product with M.
    // Computed before the solution: with DSS the solver overwrites the values of K, M stays intact.
    std :: vector< FloatArray >Mr(3);
    totalMass.resize(3);
    for ( int dir = 1; dir <= 3; dir++ ) {
        FloatArray r(neq);
        r.zero();
        for ( auto &dman : d->giveDofManagers() ) {
            for ( Dof *dof : * dman ) {
                if ( !dof->isPrimaryDof() || dof->giveDofID() != ( DofIDItem ) ( D_u + dir - 1 ) ) {
                    continue;       // slave DOFs (rigid arms) follow their master
                }
                int eq = dof->__giveEquationNumber();
                if ( eq > 0 ) {
                    r.at(eq) = 1.0;
                }
            }
        }
        M->times(r, Mr [ dir - 1 ]);
        totalMass.at(dir) = r.dotProduct(Mr [ dir - 1 ]);
    }

    auto es = classFactory.createGeneralizedEigenValueSolver(GES_LanczosSturm, d, this);
    if ( !es ) {
        OOFEM_ERROR("Lanczos eigenvalue solver not available");
    }
    n = std :: min(n, neq);
    OOFEM_LOG_INFO("ResponseSpectrum: computing %d modes\n", n);
    es->solve(* K, * M, omega2, eigVec, rtolv, n);

    // With a singular M (e.g. masses only on the DOFs of a rigid diaphragm) the modes with mass are
    // at most rank(M): only those found (omega2 > 0) are kept.
    int nv = 0;
    for ( int i = 1; i <= omega2.giveSize(); i++ ) {
        if ( omega2.at(i) > 0.0 ) {
            nv++;
        }
    }
    if ( nv < omega2.giveSize() ) {
        FloatArray o(nv);
        FloatMatrix v(neq, nv);
        for ( int i = 1, c = 0; i <= omega2.giveSize(); i++ ) {
            if ( omega2.at(i) > 0.0 ) {
                c++;
                o.at(c) = omega2.at(i);
                for ( int k = 1; k <= neq; k++ ) {
                    v.at(k, c) = eigVec.at(k, i);
                }
            }
        }
        OOFEM_LOG_INFO("ResponseSpectrum: %d modes with mass (%d requested)\n", nv, n);
        omega2 = o;
        eigVec = v;
    }
    int nm = omega2.giveSize();
    gamma.resize(nm, 3);
    for ( int i = 1; i <= nm; i++ ) {
        for ( int dir = 1; dir <= 3; dir++ ) {
            double g = 0.0;
            for ( int k = 1; k <= neq; k++ ) {
                g += eigVec.at(k, i) * Mr [ dir - 1 ].at(k);
            }
            gamma.at(i, dir) = g;
        }
    }
}

void
ResponseSpectrum :: solveYourself()
{
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    TimeStep *tStep = this->giveNextStep();
    this->updateAttributes(this->giveCurrentMetaStep() );
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering() );

    if ( pdelta ) {
        // gravity load state (loads of the input at t = 1), for the initial stress matrix
        OOFEM_LOG_INFO("ResponseSpectrum: static solution of the gravity loads for second order effects\n");
        activeMode = 0;
        uStatic.clear();
        auto K = classFactory.createSparseMtrx(sparseMtrxType);
        K->buildInternalStructure(this, 1, EModelDefaultEquationNumbering() );
        this->assemble(* K, tStep, TangentAssembler(TangentStiffness), EModelDefaultEquationNumbering(), d);
        FloatArray F(neq);
        F.zero();
        this->assembleVector(F, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), d);
        FloatArray u(neq);
        u.zero();
        auto ls = classFactory.createSparseLinSolver(solverType, d, this);
        if ( ls->solve(* K, F, u) != CR_CONVERGED ) {
            OOFEM_ERROR("static solution of the gravity loads failed");
        }
        uStatic = u;
    }

    int n = nroot;
    for ( ;; ) {
        this->computeModes(tStep, n);
        int nm = omega2.giveSize();
        bool ok = true;
        for ( int dir : directions ) {
            double s = 0.0;
            for ( int i = 1; i <= nm; i++ ) {
                s += gamma.at(i, dir) * gamma.at(i, dir);
            }
            if ( totalMass.at(dir) > 0.0 && s / totalMass.at(dir) < minMassRatio ) {
                ok = false;
            }
        }
        bool exhausted = nm < n;     // all modes with mass have already been found
        if ( ok || exhausted || n >= maxNroot || n >= neq ) {
            if ( !ok ) {
                OOFEM_WARNING("minimum effective mass ratio not reached with %d modes", n);
            }
            break;
        }
        n = std :: min(2 * n, maxNroot);
        OOFEM_LOG_INFO("ResponseSpectrum: effective mass too small, increasing to %d modes\n", n);
    }
    nroot = omega2.giveSize();

    this->terminate(tStep);
    OOFEM_LOG_INFO("ResponseSpectrum: solution time %.2f s\n", this->giveSolutionStepTime() );
}

void
ResponseSpectrum :: doStepOutput(TimeStep *tStep)
{
    this->printOutputAt(this->giveOutputStream(), tStep);
    fflush(this->giveOutputStream() );
}

void
ResponseSpectrum :: printOutputAt(FILE *file, TimeStep *tStep)
{
    Domain *dom = this->giveDomain(1);
    int nm = omega2.giveSize();
    const char *dirNames = "XYZ";
    fprintf(file, "\nOutput for time %.3e \n\n", 1.0);
    fprintf(file, "RSA combination %s, damping %.4f, modes %d%s\n",
            combination == 0 ? "CQC" : "SRSS", damping, nm, pdelta ? ", second order (K + KG of gravity loads)" : "");
    fprintf(file, "RSA total_mass %.6e %.6e %.6e\n", totalMass.at(1), totalMass.at(2), totalMass.at(3) );
    fprintf(file, "RSA mode  T[s]  omega2  Gx Gy Gz  mx my mz  (mass ratios)  cumulative\n");
    double cum [ 3 ] = { 0, 0, 0 };
    for ( int i = 1; i <= nm; i++ ) {
        double T = 2.0 * M_PI / sqrt(omega2.at(i) );
        double m [ 3 ];
        for ( int dir = 0; dir < 3; dir++ ) {
            m [ dir ] = totalMass.at(dir + 1) > 0 ? gamma.at(i, dir + 1) * gamma.at(i, dir + 1) / totalMass.at(dir + 1) : 0.0;
            cum [ dir ] += m [ dir ];
        }
        fprintf(file, "RSA mode %d %.6e %.8e %.6e %.6e %.6e %.5f %.5f %.5f %.5f %.5f %.5f\n", i, T, omega2.at(i),
                gamma.at(i, 1), gamma.at(i, 2), gamma.at(i, 3), m [ 0 ], m [ 1 ], m [ 2 ], cum [ 0 ], cum [ 1 ], cum [ 2 ]);
    }

    std :: vector< double >q(nm);
    int sg;
    for ( int dir : directions ) {
        // modal amplitude: u_i = Gamma_i Sa(T_i) / omega_i^2 phi_i
        std :: vector< double >c(nm);
        for ( int i = 1; i <= nm; i++ ) {
            c [ i - 1 ] = gamma.at(i, dir) * this->giveSpectralAcceleration(2.0 * M_PI / sqrt(omega2.at(i) ) ) / omega2.at(i);
        }
        for ( int i = 1; i <= nm; i++ ) {
            q [ i - 1 ] = c [ i - 1 ] * omega2.at(i) * gamma.at(i, dir);
        }
        double V = this->combineModalValues(q, sg);
        fprintf(file, "\nRSA direction %c\n", dirNames [ dir - 1 ]);
        fprintf(file, "RSA base_shear %.6e\n", V);

        // nodal displacements
        IntArray ids = { D_u, D_v, D_w, R_u, R_v, R_w };
        std :: vector< std :: vector< double > >val(6, std :: vector< double >(nm, 0.0) );
        for ( auto &dman : dom->giveDofManagers() ) {
            for ( int i = 1; i <= nm; i++ ) {
                activeMode = i;
                for ( int k = 0; k < 6; k++ ) {
                    val [ k ] [ i - 1 ] = dman->hasDofID( ( DofIDItem ) ids [ k ]) ?
                                          c [ i - 1 ] * dman->giveDofWithID(ids [ k ])->giveUnknown(VM_Total, tStep) : 0.0;
                }
            }
            fprintf(file, "RSA node %d", dman->giveLabel() );
            for ( int k = 0; k < 6; k++ ) {
                double v = this->combineModalValues(val [ k ], sg);
                fprintf(file, " %.10e", sg * v);
            }
            fprintf(file, "\n");
        }

        // elements
        for ( auto &el : dom->giveElements() ) {
            if ( auto *sh = dynamic_cast< DKMShell * >( el.get() ) ) {
                int nn = sh->giveNumberOfDofManagers();
                std :: vector< std :: vector< double > >v(nn * DKMSHELL_NS, std :: vector< double >(nm) );
                double S [ 4 * DKMSHELL_NS ], e [ 9 ];
                for ( int i = 1; i <= nm; i++ ) {
                    activeMode = i;
                    sh->computeVertexStresses(tStep, S);
                    for ( int k = 0; k < nn * DKMSHELL_NS; k++ ) {
                        v [ k ] [ i - 1 ] = c [ i - 1 ] * S [ k ];
                    }
                }
                sh->giveLocalAxes(e);
                fprintf(file, "RSA shell %d axes", el->giveLabel() );
                for ( int k = 0; k < 9; k++ ) {
                    fprintf(file, " %.10e", e [ k ]);
                }
                fprintf(file, "\n");
                for ( int a = 0; a < nn; a++ ) {
                    fprintf(file, "RSA vertex %d node %d", a + 1, sh->giveNode(a + 1)->giveLabel() );
                    for ( int k = 0; k < DKMSHELL_NS; k++ ) {
                        double w = this->combineModalValues(v [ a * DKMSHELL_NS + k ], sg);
                        fprintf(file, " %.10e", sg * w);
                    }
                    fprintf(file, "\n");
                }
                continue;
            }
            // other elements: end forces in the local system, f = K u
            FloatMatrix Ke;
            el->giveCharacteristicMatrix(Ke, TangentStiffnessMatrix, tStep);
            if ( Ke.giveNumberOfRows() == 0 || Ke.computeFrobeniusNorm() == 0.0 ) {
                continue;                   // e.g. lumped masses
            }
            int nd = Ke.giveNumberOfRows();
            std :: vector< std :: vector< double > >v(nd, std :: vector< double >(nm) );
            FloatArray u, f;
            for ( int i = 1; i <= nm; i++ ) {
                activeMode = i;
                el->computeVectorOf(VM_Total, tStep, u);
                f.beProductOf(Ke, u);
                for ( int k = 0; k < nd; k++ ) {
                    v [ k ] [ i - 1 ] = c [ i - 1 ] * f [ k ];
                }
            }
            fprintf(file, "RSA forces %d %d", el->giveLabel(), nd);
            for ( int k = 0; k < nd; k++ ) {
                double w = this->combineModalValues(v [ k ], sg);
                fprintf(file, " %.10e", sg * w);
            }
            fprintf(file, "\n");
        }
    }
    activeMode = 0;
    double utsec = this->timer.getUtime(EngngModelTimer :: EMTT_AnalysisTimer);
    fprintf(file, "\nUser time consumed by solution step: %.3f [s]\n\n", utsec);
}
} // end namespace oofem

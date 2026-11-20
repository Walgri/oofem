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

/*
 * Outline of the algorithm:
 *   1. K - sigma0 M is factorized (sigma0 = 0 unless OOFEM_LANCZOS_SHIFT is set) and its inertia n0 counted;
 *   2. m steps of Lanczos on the operator (K - sigma0 M)^-1 M, starting from a random vector M-orthogonal to
 *      the modes already accepted, with full reorthogonalization (up to two passes) in the M inner product;
 *   3. eigenvalues of the tridiagonal matrix (implicit QL), residual |beta_m s_mi|, acceptance of the
 *      converged pairs: lambda = sigma0 + 1 / theta;
 *   4. once there are at least nroot modes: sigma_c between the last requested and the next one; K - sigma_c M
 *      is factorized and its negative pivots N counted; if N exceeds the modes found below sigma_c the
 *      algorithm restarts from step 2 (the new vectors are orthogonal to the modes found, so they converge
 *      to the missing ones).
 */
#include "lanczossturm.h"
#include "dssmatrix.h"
#include "skyline.h"
#include "classfactory.h"
#include "error.h"
#include "logger.h"

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <memory>
#include <cstdlib>
#include <chrono>

namespace oofem {
REGISTER_GeneralizedEigenValueSolver(LanczosSturm, GES_LanczosSturm);

LanczosSturm :: LanczosSturm(Domain *d, EngngModel *m) : SparseGeneralEigenValueSystemNM(d, m) { }

namespace {

// ---------------------------------------------------------------- factorization of K - sigma M
class ShiftInvertOperator
{
public:
    SparseMtrx &K, &M;
    DSSMatrix *kd;                       // K as DSS (factorized in place, values saved)
    std :: unique_ptr< SparseMtrx > ks;  // unfactorized skyline copy of K
    std :: unique_ptr< SparseMtrx > fs;  // factorized skyline of K - sigma M
    double sigma;
    int nFactorizations;

    ShiftInvertOperator(SparseMtrx &k, SparseMtrx &m) : K(k), M(m), sigma(0.0), nFactorizations(0)
    {
        kd = dynamic_cast< DSSMatrix * >( & K );
        if ( kd ) {
            if ( !dynamic_cast< DSSMatrix * >( & M ) ) {
                OOFEM_ERROR("K and M must have the same matrix type");
            }
            kd->saveValues();
        } else if ( dynamic_cast< Skyline * >( & K ) ) {
            ks = K.clone();
        } else {
            OOFEM_ERROR("only DSS (smtype 8) or skyline (smtype 0) matrices are supported");
        }
    }

    void factorize(double s)
    {
        sigma = s;
        nFactorizations++;
        if ( kd ) {
            if ( !kd->restoreValues() ) {
                OOFEM_ERROR("restoring the values of K failed");
            }
            if ( s != 0.0 && !kd->addScaled(-s, * dynamic_cast< DSSMatrix * >( & M ) ) ) {
                OOFEM_ERROR("K and M do not have the same DSS structure");
            }
            kd->factorized();
        } else {
            fs = ks->clone();
            if ( s != 0.0 ) {
                fs->add(-s, M);
            }
            fs->factorized();
        }
    }

    /// x = (K - sigma M)^-1 b
    void solve(FloatArray &b, FloatArray &x)
    {
        if ( kd ) {
            kd->solve(b, x);
        } else {
            x = b;
            fs->backSubstitutionWith(x);
        }
    }

    /// Number of negative pivots of the current factorization
    long giveInertia()
    {
        if ( kd ) {
            return kd->countNegativePivots();
        }
        long n = 0;
        int neq = fs->giveNumberOfRows();
        for ( int i = 1; i <= neq; i++ ) {
            if ( static_cast< const SparseMtrx * >( fs.get() )->at(i, i) < 0.0 ) {
                n++;
            }
        }
        return n;
    }
};

// ---------------------------------------------------------------- implicit QL on a symmetric tridiagonal matrix
// d: diagonal (n), e: subdiagonal with e[i] = T(i, i-1), e[0] arbitrary; z (n x n, by rows) = I on input,
// on output its columns are the eigenvectors.
bool tqli(std :: vector< double > &d, std :: vector< double > &e, std :: vector< double > &z, int n)
{
    for ( int i = 1; i < n; i++ ) {
        e [ i - 1 ] = e [ i ];
    }
    e [ n - 1 ] = 0.0;
    for ( int l = 0; l < n; l++ ) {
        int iter = 0, m;
        do {
            for ( m = l; m < n - 1; m++ ) {
                double dd = fabs(d [ m ]) + fabs(d [ m + 1 ]);
                if ( fabs(e [ m ]) <= 1.e-15 * dd ) {
                    break;
                }
            }
            if ( m != l ) {
                if ( iter++ == 100 ) {
                    return false;
                }
                double g = ( d [ l + 1 ] - d [ l ] ) / ( 2.0 * e [ l ] );
                double r = hypot(g, 1.0);
                g = d [ m ] - d [ l ] + e [ l ] / ( g + ( g >= 0.0 ? fabs(r) : -fabs(r) ) );
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for ( i = m - 1; i >= l; i-- ) {
                    double f = s * e [ i ], b = c * e [ i ];
                    e [ i + 1 ] = ( r = hypot(f, g) );
                    if ( r == 0.0 ) {
                        d [ i + 1 ] -= p;
                        e [ m ] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d [ i + 1 ] - p;
                    r = ( d [ i ] - g ) * s + 2.0 * c * b;
                    d [ i + 1 ] = g + ( p = s * r );
                    g = c * r - b;
                    for ( int k = 0; k < n; k++ ) {
                        f = z [ k * n + i + 1 ];
                        z [ k * n + i + 1 ] = s * z [ k * n + i ] + c * f;
                        z [ k * n + i ] = c * z [ k * n + i ] - s * f;
                    }
                }
                if ( r == 0.0 && i >= l ) {
                    continue;
                }
                d [ l ] -= p;
                e [ l ] = g;
                e [ m ] = 0.0;
            }
        } while ( m != l );
    }
    return true;
}

// number of contiguous converged Ritz pairs starting from the largest theta (smallest lambda)
int countConvergedLeading(const std :: vector< double > &d, const std :: vector< double > &z, int m, double beta, double tol)
{
    std :: vector< int > ord(m);
    for ( int i = 0; i < m; i++ ) {
        ord [ i ] = i;
    }
    std :: sort(ord.begin(), ord.end(), [ & ](int x, int y) { return d [ x ] > d [ y ]; });
    int c = 0;
    for ( int k = 0; k < m; k++ ) {
        int i = ord [ k ];
        if ( d [ i ] <= 0.0 || fabs(beta * z [ ( m - 1 ) * m + i ]) > tol * d [ i ] ) {
            break;
        }
        c++;
    }
    return c;
}

// Symmetric matrix in compact form (one half + diagonal): product y = A x
struct CompactSymMatrix {
    std :: vector< unsigned long > cp, ri;
    std :: vector< double > v;
    bool active = false;
    void times(const FloatArray &x, FloatArray &y) const
    {
        long n = ( long ) cp.size() - 1;
        y.resize(n);
        y.zero();
        const double *px = x.givePointer();
        double *py = y.givePointer();
        for ( long j = 0; j < n; j++ ) {
            double xj = px [ j ], s = 0.0;
            for ( unsigned long k = cp [ j ]; k < cp [ j + 1 ]; k++ ) {
                long i = ( long ) ri [ k ];
                double a = v [ k ];
                if ( i == j ) {
                    s += a * xj;
                } else {
                    py [ i ] += a * xj;
                    s += a * px [ i ];
                }
            }
            py [ j ] += s;
        }
    }
};

struct EigenPair {
    double lambda;
    FloatArray y, My;
};

} // namespace


ConvergedReason
LanczosSturm :: solve(SparseMtrx &K, SparseMtrx &M, FloatArray &_eigv, FloatMatrix &_r, double rtol, int nroot)
{
    auto t0 = std :: chrono :: steady_clock :: now();
    int n = K.giveNumberOfRows();
    if ( nroot > n ) {
        nroot = n;
    }
    double sigma0 = getenv("OOFEM_LANCZOS_SHIFT") ? atof( getenv("OOFEM_LANCZOS_SHIFT") ) : 0.0;
    double tolres = std :: max(sqrt(rtol), 1.e-10);      // relative residual for accepting Ritz pairs
    int maxRounds = 20;

    // Product with M: with DSS the matrix is stored in the structure of the factor (including fill-in),
    // so a compact copy with the assembled entries only is extracted
    CompactSymMatrix Mc;
    if ( auto *md = dynamic_cast< DSSMatrix * >( & M ) ) {
        md->extractCompressed(Mc.cp, Mc.ri, Mc.v);
        Mc.active = true;
    }
    auto Mt = [ & ](const FloatArray &x, FloatArray &y) {
        if ( Mc.active ) {
            Mc.times(x, y);
        } else {
            M.times(x, y);
        }
    };
    ShiftInvertOperator op(K, M);
    op.factorize(sigma0);
    long n0 = op.giveInertia();      // eigenvalues below sigma0 (0 if K is positive definite and sigma0 = 0)
    if ( n0 > 0 ) {
        OOFEM_LOG_INFO("Lanczos: %ld eigenvalues below the shift %g (excluded from the count)\n", n0, sigma0);
    }

    std :: vector< EigenPair > accepted;
    // steps before the first convergence check; then the check is repeated every "checkInterval" steps
    int firstCheck = std :: min(n, std :: max(2 * nroot, nroot + 20));
    // maximum steps per round: limits the memory of the Lanczos vectors (m x n doubles)
    int maxSteps = std :: min( n, std :: max(6 * nroot + 60, 300) );
    if ( getenv("OOFEM_LANCZOS_MAXSTEPS") && atoi( getenv("OOFEM_LANCZOS_MAXSTEPS") ) > 0 ) {
        maxSteps = std :: min( n, atoi( getenv("OOFEM_LANCZOS_MAXSTEPS") ) );
    }
    long totalSteps = 0;
    double t_sol = 0, t_M = 0, t_ort = 0;
    long nSecondPass = 0;
    bool complete = false;
    long sturm = -1, found = 0;
    double sigmac = 0.0;
    int wantedInRound = nroot;       // converged pairs to obtain in the round

    for ( int round = 0; round < maxRounds; round++ ) {
        int mlim = std :: min(maxSteps, n - ( int ) accepted.size() );
        if ( mlim <= 0 ) {
            break;
        }
        // -------- random starting vector, filtered by the operator and M-orthogonal to the accepted modes
        std :: mt19937 gen(12345 + 7919 * round);
        std :: uniform_real_distribution< double > dist(-0.5, 0.5);
        FloatArray r(n), Mr;
        for ( int i = 1; i <= n; i++ ) {
            r.at(i) = dist(gen);
        }
        {
            FloatArray tmp;
            Mt(r, tmp);
            op.solve(tmp, r);
        }
        for ( int pass = 0; pass < 2; pass++ ) {
            Mt(r, Mr);
            for ( auto &md : accepted ) {
                r.add(-Mr.dotProduct(md.y), md.y);
            }
        }
        Mt(r, Mr);
        double beta = sqrt( std :: max(r.dotProduct(Mr), 0.0) );
        if ( beta <= 0.0 ) {
            OOFEM_WARNING("zero starting vector");
            break;
        }
        // Only the Lanczos vectors q_j are stored; products M w are recomputed (M is cheap compared
        // to a solution), so memory stays m x n instead of 2 m x n.
        std :: vector< FloatArray > Q;
        Q.reserve(mlim + 1);
        r.times(1.0 / beta);
        Mr.times(1.0 / beta);
        Q.push_back(r);
        FloatArray Mq = Mr;              // M q_j of the current vector
        std :: vector< double > alpha, bet, d, z;
        int meff = 0, nextCheck = std :: min(firstCheck, mlim), checkInterval = std :: max(10, nroot / 2);
        FloatArray w, Mw;
        bool done = false;
        for ( int j = 0; !done; j++ ) {
            auto c0 = std :: chrono :: steady_clock :: now();
            op.solve(Mq, w);
            auto c1 = std :: chrono :: steady_clock :: now();
            t_sol += std :: chrono :: duration< double >(c1 - c0).count();
            totalSteps++;
            if ( j > 0 ) {
                w.add(-bet [ j - 1 ], Q [ j - 1 ]);
            }
            double a = w.dotProduct(Mq);
            alpha.push_back(a);
            w.add(-a, Q [ j ]);
            // full reorthogonalization in M (classical Gram-Schmidt); second pass only if the first
            // one reduced the norm considerably (DGKS criterion, threshold 1/sqrt(2))
            Mt(w, Mw);
            for ( int pass = 0; pass < 2; pass++ ) {
                double n2before = w.dotProduct(Mw);
                auto c3 = std :: chrono :: steady_clock :: now();
                {
                    double *pw = w.givePointer();
                    const double *pm = Mw.givePointer();
                    for ( int i = 0; i <= j; i++ ) {
                        const double *pq = Q [ i ].givePointer();
                        double c = 0.0;
                        for ( int k = 0; k < n; k++ ) {
                            c += pm [ k ] * pq [ k ];
                        }
                        for ( int k = 0; k < n; k++ ) {
                            pw [ k ] -= c * pq [ k ];
                        }
                    }
                }
                for ( auto &md : accepted ) {
                    w.add(-Mw.dotProduct(md.y), md.y);
                }
                auto c4 = std :: chrono :: steady_clock :: now();
                t_ort += std :: chrono :: duration< double >(c4 - c3).count();
                Mt(w, Mw);
                t_M += std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - c4).count();
                if ( w.dotProduct(Mw) > 0.5 * n2before ) {
                    break;
                }
                nSecondPass++;
            }
            double b = sqrt( std :: max(w.dotProduct(Mw), 0.0) );
            bet.push_back(b);
            meff = j + 1;
            bool breakdown = b <= 1.e-12 * fabs(a);
            if ( meff >= nextCheck || meff >= mlim || breakdown ) {
                // -------- eigenvalues of the tridiagonal matrix and count of the converged pairs
                d.assign(alpha.begin(), alpha.end());
                std :: vector< double > e(meff, 0.0);
                for ( int i = 1; i < meff; i++ ) {
                    e [ i ] = bet [ i - 1 ];
                }
                z.assign(meff * meff, 0.0);
                for ( int i = 0; i < meff; i++ ) {
                    z [ i * meff + i ] = 1.0;
                }
                if ( !tqli(d, e, z, meff) ) {
                    OOFEM_WARNING("QL iteration on the tridiagonal matrix did not converge");
                }
                // leading converged pairs: decreasing theta (increasing lambda), contiguous from the first;
                // pairs at the other end of the spectrum converge early but are of no interest
                int conv = countConvergedLeading(d, z, meff, b, tolres);
                if ( conv >= wantedInRound || meff >= mlim || breakdown ) {
                    done = true;
                } else {
                    nextCheck = std :: min(mlim, meff + checkInterval);
                }
            }
            if ( !done ) {
                w.times(1.0 / b);
                Mw.times(1.0 / b);
                Q.push_back(w);
                Mq = Mw;
            }
        }
        // -------- acceptance of the converged Ritz pairs
        double betam = bet [ meff - 1 ];
        int nNew = 0;
        std :: vector< int > ord(meff);
        for ( int i = 0; i < meff; i++ ) {
            ord [ i ] = i;
        }
        std :: sort(ord.begin(), ord.end(), [ & ](int x, int y) { return d [ x ] > d [ y ]; });
        int nLeading = countConvergedLeading(d, z, meff, betam, tolres);
        for ( int k = 0; k < nLeading; k++ ) {
            int i = ord [ k ];
            double theta = d [ i ];
            EigenPair md;
            md.lambda = sigma0 + 1.0 / theta;
            md.y.resize(n);
            md.y.zero();
            for ( int j = 0; j < meff; j++ ) {
                md.y.add(z [ j * meff + i ], Q [ j ]);
            }
            // purification (Nour-Omid, Parlett, Ericsson, Jensen 1987): with a singular M the massless DOFs
            // are not controlled by the M norm and accumulate round-off; one application of the operator,
            // y <- (K - sigma M)^-1 M y / theta, rebuilds them by static condensation
            {
                FloatArray My, yp;
                Mt(md.y, My);
                op.solve(My, yp);
                yp.times(1.0 / theta);           // (K - sigma M)^-1 M y = theta y
                md.y = yp;
                totalSteps++;
            }
            // reorthogonalize against the accepted modes and normalize in M
            Mt(md.y, md.My);
            for ( auto &ac : accepted ) {
                md.y.add(-md.My.dotProduct(ac.y), ac.y);
            }
            Mt(md.y, md.My);
            double nm = sqrt( std :: max(md.y.dotProduct(md.My), 0.0) );
            if ( nm < 0.5 ) {
                continue;   // vector already contained in the accepted modes
            }
            md.y.times(1.0 / nm);
            md.My.clear();  // not needed any more
            accepted.push_back(std :: move(md) );
            nNew++;
        }
        Q.clear();
        std :: sort(accepted.begin(), accepted.end(), [](const EigenPair &a, const EigenPair &b) { return a.lambda < b.lambda; });
        OOFEM_LOG_INFO("Lanczos: round %d, %d steps, %d new pairs accepted, total %d\n",
                       round + 1, meff, nNew, ( int ) accepted.size() );

        if ( ( int ) accepted.size() < nroot ) {
            wantedInRound = nroot - ( int ) accepted.size();
            continue;
        }
        // -------- Sturm sequence check
        double lp = accepted [ nroot - 1 ].lambda;
        sigmac = lp * ( 1.0 + 1.e-6 ) + 1.e-12;
        if ( ( int ) accepted.size() > nroot ) {
            double lq = accepted [ nroot ].lambda;
            if ( lq > lp * ( 1.0 + 1.e-6 ) ) {
                sigmac = 0.5 * ( lp + lq );
            }
        }
        op.factorize(sigmac);
        sturm = op.giveInertia() - n0;
        found = 0;
        for ( auto &md : accepted ) {
            if ( md.lambda < sigmac ) {
                found++;
            }
        }
        OOFEM_LOG_INFO("Sturm check: sigma = %.10e, eigenvalues below sigma = %ld, found = %ld, requested = %d -> %s\n",
                       sigmac, sturm, found, nroot, sturm == found ? "complete" : "MISSING MODES");
        if ( sturm == found ) {
            complete = true;
            break;
        }
        // modes missing below sigma_c (typically multiple eigenvalues): new round orthogonal to those found
        op.factorize(sigma0);
        wantedInRound = ( int ) std :: min( ( long ) nroot, std :: max(1L, sturm - found) );
    }

    OOFEM_LOG_INFO("Lanczos: times: solutions %.1f s, M products (reorth.) %.1f s, reorthogonalization %.1f s (second pass in %ld steps)\n", t_sol, t_M, t_ort, nSecondPass);
    double dt = std :: chrono :: duration< double >(std :: chrono :: steady_clock :: now() - t0).count();
    OOFEM_LOG_INFO("Lanczos: %s; %ld Lanczos steps, %d factorizations, %.1f s\n",
                   complete ? "set of modes complete (Sturm check)" : "WARNING: completeness NOT verified",
                   totalSteps, op.nFactorizations, dt);

    int nout = std :: min(nroot, ( int ) accepted.size() );
    _eigv.resize(nroot);
    _eigv.zero();
    _r.resize(n, nroot);
    _r.zero();
    for ( int i = 1; i <= nout; i++ ) {
        _eigv.at(i) = accepted [ i - 1 ].lambda;
        for ( int j = 1; j <= n; j++ ) {
            _r.at(j, i) = accepted [ i - 1 ].y.at(j);
        }
    }
    return complete ? CR_CONVERGED : CR_DIVERGED_ITS;
}
} // end namespace oofem

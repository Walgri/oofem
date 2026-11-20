/*
 * Computational kernel of the flat shell elements DKMTShell and DKMQShell: ITW membrane +
 * DKMT/DKMQ plate. See dkmshellkernel.h.
 *
 * References:
 *   D.J. Allman, "A compatible triangular element including vertex rotations for plane
 *     elasticity analysis", Comput. Struct. 19 (1984).
 *   A. Ibrahimbegovic, R.L. Taylor, E.L. Wilson, "A robust quadrilateral membrane finite
 *     element with drilling degrees of freedom", IJNME 30 (1990) 445-457.
 *   I. Katili, "A new discrete Kirchhoff-Mindlin element based on Mindlin-Reissner plate
 *     theory and assumed shear strain fields", IJNME 36 (1993), part I (DKMT) and II (DKMQ).
 *
 * Written to be read: explicit loops, no libraries.
 */
#include "dkmshellkernel.h"
#include <math.h>
#include <string.h>

#define NMAX 4
#define DMAX 24

/* local DOFs per node */
#define GU 0
#define GV 1
#define GW 2
#define GRX 3
#define GRY 4
#define GRZ 5

/* ------------------------------------------------------------------ element data */
typedef struct {
    int nn, nd;
    double x[NMAX], y[NMAX];
    double L[NMAX], C[NMAX], S[NMAX];     /* side k: from node k to node k+1; length, direction cosines */
    double phi[NMAX];                     /* Katili shear parameter per side */
    double A[NMAX][DMAX];                 /* delta-beta_k (quadratic tangential rotation at mid-side) */
    double G[NMAX][DMAX];                 /* gamma_s,k (tangential shear strain, constant along the side) */
    double Rinv[3][3];                    /* triangle: linear shear field (Nedelec) from the three gamma_s */
    double xc, yc;
    double E, nu, t, k, Gm, gdrill;
} ElementData;

/* integration / evaluation point */
typedef struct {
    double N[NMAX];
    double detJ, Jinv[2][2];
    double xp, yp;
    double Bm[3][DMAX];     /* membrane strains: ex ey gxy */
    double bsk[DMAX];       /* skew-symmetric part of the gradient minus the rotation: 1/2(v,x-u,y) - omega */
    double Bb[3][DMAX];     /* curvatures: bx,x  by,y  bx,y+by,x */
    double Bs[2][DMAX];     /* transverse shear strains: gxz gyz (assumed field) */
    double gu[2][DMAX], gv[2][DMAX], gw[2][DMAX], gwl[2][DMAX];   /* gradients for the initial stress matrix */
} PointData;

static void shape_functions(int nn, double xi, double eta, double *N, double *dNxi, double *dNeta,
                     double *P, double *dPxi, double *dPeta)
{
    if ( nn == 3 ) {
        double l1 = 1.0 - xi - eta, l2 = xi, l3 = eta;
        N[0] = l1; N[1] = l2; N[2] = l3;
        dNxi[0] = -1.0; dNeta[0] = -1.0;
        dNxi[1] =  1.0; dNeta[1] =  0.0;
        dNxi[2] =  0.0; dNeta[2] =  1.0;
        /* side bubbles: value 1 at mid-side */
        P[0] = 4.0 * l1 * l2;  dPxi[0] = 4.0 * ( 1.0 - 2.0 * xi - eta );  dPeta[0] = -4.0 * xi;
        P[1] = 4.0 * l2 * l3;  dPxi[1] = 4.0 * eta;                        dPeta[1] = 4.0 * xi;
        P[2] = 4.0 * l3 * l1;  dPxi[2] = -4.0 * eta;                       dPeta[2] = 4.0 * ( 1.0 - xi - 2.0 * eta );
    } else {
        static const double xn[4] = { -1.0, 1.0, 1.0, -1.0 }, en[4] = { -1.0, -1.0, 1.0, 1.0 };
        for ( int i = 0; i < 4; i++ ) {
            N[i] = 0.25 * ( 1.0 + xi * xn[i] ) * ( 1.0 + eta * en[i] );
            dNxi[i] = 0.25 * xn[i] * ( 1.0 + eta * en[i] );
            dNeta[i] = 0.25 * en[i] * ( 1.0 + xi * xn[i] );
        }
        P[0] = 0.5 * ( 1.0 - xi * xi ) * ( 1.0 - eta );   dPxi[0] = -xi * ( 1.0 - eta );            dPeta[0] = -0.5 * ( 1.0 - xi * xi );
        P[1] = 0.5 * ( 1.0 + xi ) * ( 1.0 - eta * eta );  dPxi[1] = 0.5 * ( 1.0 - eta * eta );      dPeta[1] = -eta * ( 1.0 + xi );
        P[2] = 0.5 * ( 1.0 - xi * xi ) * ( 1.0 + eta );   dPxi[2] = -xi * ( 1.0 + eta );            dPeta[2] = 0.5 * ( 1.0 - xi * xi );
        P[3] = 0.5 * ( 1.0 - xi ) * ( 1.0 - eta * eta );  dPxi[3] = -0.5 * ( 1.0 - eta * eta );     dPeta[3] = -eta * ( 1.0 - xi );
    }
}

static int prepare(ElementData *e, int nn, const double *x, const double *y, const DKMShellMaterial *m)
{
    memset(e, 0, sizeof( ElementData ) );
    if ( nn != 3 && nn != 4 ) {
        return -1;
    }
    e->nn = nn;
    e->nd = 6 * nn;
    e->E = m->E;
    e->nu = m->nu;
    e->t = m->t;
    e->k = m->k > 0.0 ? m->k : 5.0 / 6.0;
    e->Gm = m->E / ( 2.0 * ( 1.0 + m->nu ) );
    e->gdrill = m->gdrill > 0.0 ? m->gdrill : e->Gm;
    for ( int i = 0; i < nn; i++ ) {
        e->x[i] = x[i];
        e->y[i] = y[i];
        e->xc += x[i] / nn;
        e->yc += y[i] / nn;
    }
    for ( int k = 0; k < nn; k++ ) {
        int i = k, j = ( k + 1 ) % nn;
        double dx = x[j] - x[i], dy = y[j] - y[i];
        double L = sqrt(dx * dx + dy * dy);
        if ( L <= 0.0 ) {
            return -1;
        }
        e->L[k] = L;
        e->C[k] = dx / L;
        e->S[k] = dy / L;
        /* Katili: phi_k = 2/(k(1-nu)) (t/L_k)^2 ;  phi -> 0 gives back DKT/DKQ */
        e->phi[k] = 2.0 / ( e->k * ( 1.0 - e->nu ) ) * ( e->t / L ) * ( e->t / L );
        /* shear constraint integrated along the side (Timoshenko along the side):
         *   w_j - w_i + L/2 (bs_i + bs_j) + 2/3 L dB = L gamma_s,   gamma_s = -2/3 phi dB
         *   => dB = -3 / (2 L (1+phi)) [ w_j - w_i + L/2 (bs_i + bs_j) ]
         * with bs = C bx + S by, bx = ry, by = -rx */
        double c = -3.0 / ( 2.0 * L * ( 1.0 + e->phi[k] ) );
        double *A = e->A[k];
        A[6 * j + GW] += c;
        A[6 * i + GW] -= c;
        double h = c * 0.5 * L;
        A[6 * i + GRY] += h * e->C[k];
        A[6 * j + GRY] += h * e->C[k];
        A[6 * i + GRX] -= h * e->S[k];
        A[6 * j + GRX] -= h * e->S[k];
        for ( int d = 0; d < e->nd; d++ ) {
            e->G[k][d] = -2.0 / 3.0 * e->phi[k] * A[d];
        }
    }
    if ( nn == 3 ) {
        /* shear field g = (a1 - b (y-yc), a2 + b (x-xc)): constant tangential component on every side */
        double R[3][3];
        for ( int k = 0; k < 3; k++ ) {
            int i = k, j = ( k + 1 ) % 3;
            double mx = 0.5 * ( x[i] + x[j] ) - e->xc, my = 0.5 * ( y[i] + y[j] ) - e->yc;
            R[k][0] = e->C[k];
            R[k][1] = e->S[k];
            R[k][2] = -e->C[k] * my + e->S[k] * mx;
        }
        double det = R[0][0] * ( R[1][1] * R[2][2] - R[1][2] * R[2][1] )
                   - R[0][1] * ( R[1][0] * R[2][2] - R[1][2] * R[2][0] )
                   + R[0][2] * ( R[1][0] * R[2][1] - R[1][1] * R[2][0] );
        if ( fabs(det) < 1e-14 ) {
            return -1;
        }
        e->Rinv[0][0] = ( R[1][1] * R[2][2] - R[1][2] * R[2][1] ) / det;
        e->Rinv[0][1] = ( R[0][2] * R[2][1] - R[0][1] * R[2][2] ) / det;
        e->Rinv[0][2] = ( R[0][1] * R[1][2] - R[0][2] * R[1][1] ) / det;
        e->Rinv[1][0] = ( R[1][2] * R[2][0] - R[1][0] * R[2][2] ) / det;
        e->Rinv[1][1] = ( R[0][0] * R[2][2] - R[0][2] * R[2][0] ) / det;
        e->Rinv[1][2] = ( R[0][2] * R[1][0] - R[0][0] * R[1][2] ) / det;
        e->Rinv[2][0] = ( R[1][0] * R[2][1] - R[1][1] * R[2][0] ) / det;
        e->Rinv[2][1] = ( R[0][1] * R[2][0] - R[0][0] * R[2][1] ) / det;
        e->Rinv[2][2] = ( R[0][0] * R[1][1] - R[0][1] * R[1][0] ) / det;
    }
    return 0;
}

/* All fields at a point (xi, eta) */
static int evaluate_point(const ElementData *e, double xi, double eta, PointData *p)
{
    int nn = e->nn, nd = e->nd;
    double dNxi[NMAX], dNeta[NMAX], P[NMAX], dPxi[NMAX], dPeta[NMAX];
    memset(p, 0, sizeof( PointData ) );
    shape_functions(nn, xi, eta, p->N, dNxi, dNeta, P, dPxi, dPeta);

    double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
    for ( int i = 0; i < nn; i++ ) {
        J11 += dNxi[i] * e->x[i];
        J12 += dNxi[i] * e->y[i];
        J21 += dNeta[i] * e->x[i];
        J22 += dNeta[i] * e->y[i];
        p->xp += p->N[i] * e->x[i];
        p->yp += p->N[i] * e->y[i];
    }
    p->detJ = J11 * J22 - J12 * J21;
    if ( p->detJ <= 0.0 ) {
        return -1;          /* clockwise nodes or distorted element */
    }
    p->Jinv[0][0] = J22 / p->detJ;
    p->Jinv[0][1] = -J12 / p->detJ;
    p->Jinv[1][0] = -J21 / p->detJ;
    p->Jinv[1][1] = J11 / p->detJ;

    double dNx[NMAX], dNy[NMAX], dPx[NMAX], dPy[NMAX];
    for ( int i = 0; i < nn; i++ ) {
        dNx[i] = p->Jinv[0][0] * dNxi[i] + p->Jinv[0][1] * dNeta[i];
        dNy[i] = p->Jinv[1][0] * dNxi[i] + p->Jinv[1][1] * dNeta[i];
        dPx[i] = p->Jinv[0][0] * dPxi[i] + p->Jinv[0][1] * dPeta[i];
        dPy[i] = p->Jinv[1][0] * dPxi[i] + p->Jinv[1][1] * dPeta[i];
    }

    /* ---------------- membrane (Allman / ITW)
     *   u = sum N_i u_i + sum_k P_k (y_j - y_i)/8 (om_j - om_i)
     *   v = sum N_i v_i - sum_k P_k (x_j - x_i)/8 (om_j - om_i)
     *   om = sum N_i om_i  (independent rotation)                        */
    double ux[DMAX] = { 0 }, uy[DMAX] = { 0 }, vx[DMAX] = { 0 }, vy[DMAX] = { 0 };
    for ( int i = 0; i < nn; i++ ) {
        ux[6 * i + GU] += dNx[i];
        uy[6 * i + GU] += dNy[i];
        vx[6 * i + GV] += dNx[i];
        vy[6 * i + GV] += dNy[i];
    }
    for ( int k = 0; k < nn; k++ ) {
        int i = k, j = ( k + 1 ) % nn;
        double ax = ( e->y[j] - e->y[i] ) / 8.0, ay = -( e->x[j] - e->x[i] ) / 8.0;
        ux[6 * j + GRZ] += dPx[k] * ax;  ux[6 * i + GRZ] -= dPx[k] * ax;
        uy[6 * j + GRZ] += dPy[k] * ax;  uy[6 * i + GRZ] -= dPy[k] * ax;
        vx[6 * j + GRZ] += dPx[k] * ay;  vx[6 * i + GRZ] -= dPx[k] * ay;
        vy[6 * j + GRZ] += dPy[k] * ay;  vy[6 * i + GRZ] -= dPy[k] * ay;
    }
    for ( int d = 0; d < nd; d++ ) {
        p->Bm[0][d] = ux[d];
        p->Bm[1][d] = vy[d];
        p->Bm[2][d] = uy[d] + vx[d];
        p->bsk[d] = 0.5 * ( vx[d] - uy[d] );
        p->gu[0][d] = ux[d];
        p->gu[1][d] = uy[d];
        p->gv[0][d] = vx[d];
        p->gv[1][d] = vy[d];
    }
    for ( int i = 0; i < nn; i++ ) {
        p->bsk[6 * i + GRZ] -= p->N[i];
    }

    /* ---------------- plate (DKMT / DKMQ)
     *   bx = sum N_i bx_i + sum_k P_k C_k dB_k,   by = sum N_i by_i + sum_k P_k S_k dB_k
     *   bx = ry, by = -rx                                                   */
    double bx[DMAX] = { 0 }, by[DMAX] = { 0 }, bxx[DMAX] = { 0 }, bxy[DMAX] = { 0 }, byx[DMAX] = { 0 }, byy[DMAX] = { 0 };
    for ( int i = 0; i < nn; i++ ) {
        bx[6 * i + GRY] += p->N[i];
        bxx[6 * i + GRY] += dNx[i];
        bxy[6 * i + GRY] += dNy[i];
        by[6 * i + GRX] -= p->N[i];
        byx[6 * i + GRX] -= dNx[i];
        byy[6 * i + GRX] -= dNy[i];
    }
    for ( int k = 0; k < nn; k++ ) {
        for ( int d = 0; d < nd; d++ ) {
            double a = e->A[k][d];
            if ( a == 0.0 ) {
                continue;
            }
            bx[d] += P[k] * e->C[k] * a;
            bxx[d] += dPx[k] * e->C[k] * a;
            bxy[d] += dPy[k] * e->C[k] * a;
            by[d] += P[k] * e->S[k] * a;
            byx[d] += dPx[k] * e->S[k] * a;
            byy[d] += dPy[k] * e->S[k] * a;
        }
    }
    for ( int d = 0; d < nd; d++ ) {
        p->Bb[0][d] = bxx[d];
        p->Bb[1][d] = byy[d];
        p->Bb[2][d] = bxy[d] + byx[d];
        /* discrete Kirchhoff: w,x = -bx, w,y = -by */
        p->gw[0][d] = -bx[d];
        p->gw[1][d] = -by[d];
    }
    for ( int i = 0; i < nn; i++ ) {
        p->gwl[0][6 * i + GW] = dNx[i];
        p->gwl[1][6 * i + GW] = dNy[i];
    }

    /* ---------------- assumed shear from the gamma_s constant along the sides */
    if ( nn == 3 ) {
        double dx = p->xp - e->xc, dy = p->yp - e->yc;
        for ( int k = 0; k < 3; k++ ) {
            double cx = e->Rinv[0][k] - dy * e->Rinv[2][k];
            double cy = e->Rinv[1][k] + dx * e->Rinv[2][k];
            for ( int d = 0; d < nd; d++ ) {
                p->Bs[0][d] += cx * e->G[k][d];
                p->Bs[1][d] += cy * e->G[k][d];
            }
        }
    } else {
        /* covariant components (MITC4-like): sides 0 and 2 along xi, sides 1 and 3 along eta;
         * sides 2 (3->4) and 3 (4->1) run opposite to the coordinate */
        double cxi0 = 0.5 * ( 1.0 - eta ) * 0.5 * e->L[0];
        double cxi2 = -0.5 * ( 1.0 + eta ) * 0.5 * e->L[2];
        double cet1 = 0.5 * ( 1.0 + xi ) * 0.5 * e->L[1];
        double cet3 = -0.5 * ( 1.0 - xi ) * 0.5 * e->L[3];
        for ( int d = 0; d < nd; d++ ) {
            double gxi = cxi0 * e->G[0][d] + cxi2 * e->G[2][d];
            double get = cet1 * e->G[1][d] + cet3 * e->G[3][d];
            p->Bs[0][d] = p->Jinv[0][0] * gxi + p->Jinv[0][1] * get;
            p->Bs[1][d] = p->Jinv[1][0] * gxi + p->Jinv[1][1] * get;
        }
    }
    return 0;
}

/* ------------------------------------------------------------------ integration */
typedef struct {
    int n;
    double xi[9], eta[9], w[9];
} Rule;

static void rule_full(int nn, Rule *r)
{
    if ( nn == 3 ) {
        r->n = 3;
        r->xi[0] = 1.0 / 6.0; r->eta[0] = 1.0 / 6.0;
        r->xi[1] = 2.0 / 3.0; r->eta[1] = 1.0 / 6.0;
        r->xi[2] = 1.0 / 6.0; r->eta[2] = 2.0 / 3.0;
        r->w[0] = r->w[1] = r->w[2] = 1.0 / 6.0;
    } else {
        double g = 1.0 / sqrt(3.0);
        r->n = 4;
        r->xi[0] = -g; r->eta[0] = -g;
        r->xi[1] = g;  r->eta[1] = -g;
        r->xi[2] = g;  r->eta[2] = g;
        r->xi[3] = -g; r->eta[3] = g;
        r->w[0] = r->w[1] = r->w[2] = r->w[3] = 1.0;
    }
}

/* quadrilateral 3x3: for the membrane, whose (Allman) field is 8-node serendipity */
static void rule_3x3(Rule *r)
{
    static const double g[3] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
    static const double w[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    r->n = 9;
    for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
            r->xi[3 * i + j] = g[i];
            r->eta[3 * i + j] = g[j];
            r->w[3 * i + j] = w[i] * w[j];
        }
    }
}

static void rule_center(int nn, Rule *r)
{
    r->n = 1;
    if ( nn == 3 ) {
        r->xi[0] = r->eta[0] = 1.0 / 3.0;
        r->w[0] = 0.5;
    } else {
        r->xi[0] = r->eta[0] = 0.0;
        r->w[0] = 4.0;
    }
}

/* K += f * B^T D B (D symmetric nr x nr) */
static void add_BtDB(int nd, int nr, double B[][DMAX], const double D[3][3], double f, double *K)
{
    double DB[3][DMAX];
    for ( int a = 0; a < nr; a++ ) {
        for ( int d = 0; d < nd; d++ ) {
            double s = 0.0;
            for ( int b = 0; b < nr; b++ ) {
                s += D[a][b] * B[b][d];
            }
            DB[a][d] = s;
        }
    }
    for ( int i = 0; i < nd; i++ ) {
        for ( int j = 0; j < nd; j++ ) {
            double s = 0.0;
            for ( int a = 0; a < nr; a++ ) {
                s += B[a][i] * DB[a][j];
            }
            K[i * nd + j] += f * s;
        }
    }
}

static void plane_matrix(double E, double nu, double factor, double D[3][3])
{
    double c = factor * E / ( 1.0 - nu * nu );
    D[0][0] = c;      D[0][1] = c * nu; D[0][2] = 0.0;
    D[1][0] = c * nu; D[1][1] = c;      D[1][2] = 0.0;
    D[2][0] = 0.0;    D[2][1] = 0.0;    D[2][2] = c * 0.5 * ( 1.0 - nu );
}

double dkmshell_area(int nn, const double *x, const double *y)
{
    double a = 0.0;
    for ( int i = 0; i < nn; i++ ) {
        int j = ( i + 1 ) % nn;
        a += x[i] * y[j] - x[j] * y[i];
    }
    return 0.5 * a;
}

int dkmshell_stiffness(int nn, const double *x, const double *y, const DKMShellMaterial *m, double *K)
{
    ElementData e;
    if ( prepare(&e, nn, x, y, m) ) {
        return -1;
    }
    int nd = e.nd;
    memset(K, 0, sizeof( double ) * nd * nd);
    double Dm[3][3], Db[3][3], Ds[3][3] = { { 0 } };
    plane_matrix(e.E, e.nu, e.t, Dm);
    plane_matrix(e.E, e.nu, e.t * e.t * e.t / 12.0, Db);
    Ds[0][0] = Ds[1][1] = e.k * e.Gm * e.t;

    Rule r;
    PointData p;
    /* plate: 2x2 / 3 points */
    rule_full(nn, &r);
    for ( int g = 0; g < r.n; g++ ) {
        if ( evaluate_point(&e, r.xi[g], r.eta[g], &p) ) {
            return -1;
        }
        double dA = p.detJ * r.w[g];
        add_BtDB(nd, 3, p.Bb, Db, dA, K);
        add_BtDB(nd, 2, p.Bs, Ds, dA, K);
    }
    /* membrane: triangle 3 points (exact); quadrilateral 3x3, because the Allman field is serendipity
     * and 2x2 would leave the spurious mode of the reduced-integrated Q8 */
    if ( nn == 4 ) {
        rule_3x3(&r);
    }
    for ( int g = 0; g < r.n; g++ ) {
        if ( evaluate_point(&e, r.xi[g], r.eta[g], &p) ) {
            return -1;
        }
        add_BtDB(nd, 3, p.Bm, Dm, p.detJ * r.w[g], K);
    }
    /* drilling penalty (ITW: gamma = G). Triangle: one point, as Allman/Ibrahimbegovic;
     * quadrilateral: 2x2, because one point leaves an hourglass mode on rectangles. */
#ifdef DKMSHELL_DRILL_ONE_POINT
    rule_center(nn, &r);
#else
    if ( nn == 3 ) {
        rule_center(nn, &r);
    } else {
        rule_full(nn, &r);
    }
#endif
    for ( int g = 0; g < r.n; g++ ) {
        if ( evaluate_point(&e, r.xi[g], r.eta[g], &p) ) {
            return -1;
        }
        double f = e.gdrill * e.t * p.detJ * r.w[g];
        for ( int i = 0; i < nd; i++ ) {
            for ( int j = 0; j < nd; j++ ) {
                K[i * nd + j] += f * p.bsk[i] * p.bsk[j];
            }
        }
    }
    return 0;
}

static double dot_product(int nd, const double *a, const double *u)
{
    double s = 0.0;
    for ( int d = 0; d < nd; d++ ) {
        s += a[d] * u[d];
    }
    return s;
}

int dkmshell_initial_stress(int nn, const double *x, const double *y, const DKMShellMaterial *m,
                            const double *ul, int option, double *KG)
{
    ElementData e;
    if ( prepare(&e, nn, x, y, m) ) {
        return -1;
    }
    int nd = e.nd;
    memset(KG, 0, sizeof( double ) * nd * nd);
    double Dm[3][3];
    plane_matrix(e.E, e.nu, e.t, Dm);
    Rule r;
    PointData p;
    rule_full(nn, &r);
    for ( int g = 0; g < r.n; g++ ) {
        if ( evaluate_point(&e, r.xi[g], r.eta[g], &p) ) {
            return -1;
        }
        double eps[3], Nf[3];
        for ( int a = 0; a < 3; a++ ) {
            eps[a] = dot_product(nd, p.Bm[a], ul);
        }
        for ( int a = 0; a < 3; a++ ) {
            Nf[a] = Dm[a][0] * eps[0] + Dm[a][1] * eps[1] + Dm[a][2] * eps[2];
        }
        double Sg[3][3] = { { Nf[0], Nf[2], 0 }, { Nf[2], Nf[1], 0 }, { 0, 0, 0 } };
        double dA = p.detJ * r.w[g];
        add_BtDB(nd, 2, p.gu, Sg, dA, KG);
        add_BtDB(nd, 2, p.gv, Sg, dA, KG);
        if ( option == DKMSHELL_KG_W_LINEAR ) {
            add_BtDB(nd, 2, p.gwl, Sg, dA, KG);
        } else {
            add_BtDB(nd, 2, p.gw, Sg, dA, KG);
        }
    }
    return 0;
}

int dkmshell_point_stresses(int nn, const double *x, const double *y, const DKMShellMaterial *m,
                            const double *ul, double xi, double eta, double *S)
{
    ElementData e;
    if ( prepare(&e, nn, x, y, m) ) {
        return -1;
    }
    int nd = e.nd;
    double Dm[3][3], Db[3][3];
    plane_matrix(e.E, e.nu, e.t, Dm);
    plane_matrix(e.E, e.nu, e.t * e.t * e.t / 12.0, Db);
    PointData p;
    if ( evaluate_point(&e, xi, eta, &p) ) {
        return -1;
    }
    double eps[3], kap[3];
    for ( int a = 0; a < 3; a++ ) {
        eps[a] = dot_product(nd, p.Bm[a], ul);
        kap[a] = dot_product(nd, p.Bb[a], ul);
    }
    for ( int a = 0; a < 3; a++ ) {
        S[a] = Dm[a][0] * eps[0] + Dm[a][1] * eps[1] + Dm[a][2] * eps[2];
        S[3 + a] = Db[a][0] * kap[0] + Db[a][1] * kap[1] + Db[a][2] * kap[2];
    }
    double ds = e.k * e.Gm * e.t;
    S[6] = ds * dot_product(nd, p.Bs[0], ul);
    S[7] = ds * dot_product(nd, p.Bs[1], ul);

    return 0;
}

int dkmshell_vertex_stresses(int nn, const double *x, const double *y, const DKMShellMaterial *m,
                             const double *ul, double *S)
{
    static const double vt[3][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 } };
    static const double vq[4][2] = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };
    for ( int i = 0; i < nn; i++ ) {
        double xi = nn == 3 ? vt[i][0] : vq[i][0];
        double eta = nn == 3 ? vt[i][1] : vq[i][1];
        if ( dkmshell_point_stresses(nn, x, y, m, ul, xi, eta, S + DKMSHELL_NS * i) ) {
            return -1;
        }
    }
    return 0;
}

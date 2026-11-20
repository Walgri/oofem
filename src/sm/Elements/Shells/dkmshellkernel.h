/*
 * Computational kernel of the flat shell elements DKMTShell and DKMQShell (3 or 4 nodes, 6 DOFs per node).
 *
 *   membrane: Allman / Ibrahimbegovic-Taylor-Wilson (ITW), true drilling rotation with a penalty
 *             on the skew-symmetric part of the displacement gradient (gamma = G);
 *   plate:    DKMT (3 nodes) / DKMQ (4 nodes) by Katili, discrete Kirchhoff-Mindlin:
 *             free of shear locking, transverse shear obtained directly.
 *
 * Everything is in the element local system: nodes in the x-y plane, numbered counterclockwise
 * when seen from the local z axis. DOFs per node, in this order:
 *     u  v  w  rx  ry  rz        (rx, ry, rz: rotations about the local axes, right-hand rule)
 * Matrices are square of order 6*nn, stored by rows.
 *
 * Plain C with no dependencies other than libm, so the kernel can also be built as a shared
 * library and called outside OOFEM (e.g. from Python through ctypes).
 */
#ifndef dkmshellkernel_h
#define dkmshellkernel_h

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double E;          /* Young's modulus */
    double nu;         /* Poisson's ratio */
    double t;          /* thickness */
    double k;          /* shear correction factor (5/6) */
    double gdrill;     /* drilling penalty modulus; <= 0 -> G */
} DKMShellMaterial;

/* Options of the initial stress matrix */
#define DKMSHELL_KG_W_ROTATIONS 0   /* gradient of w taken from the rotations (DK beta field), default */
#define DKMSHELL_KG_W_LINEAR    1   /* gradient of w from the nodal w only (linear/bilinear functions) */

/* Stiffness matrix K (6nn x 6nn). Returns 0 on success, -1 if the geometry is degenerate. */
int dkmshell_stiffness(int nn, const double *x, const double *y, const DKMShellMaterial *m, double *K);

/* Initial stress matrix KG due to the membrane forces produced by the local displacements ul (6nn).
 * Critical load: (K + lambda KG) phi = 0. */
int dkmshell_initial_stress(int nn, const double *x, const double *y, const DKMShellMaterial *m,
                            const double *ul, int option, double *KG);

/* Number of stress resultant components per point */
#define DKMSHELL_NS 8

/* Stress resultants at the vertices from the element field, per unit width, in the local system.
 * S[DKMSHELL_NS*i + c], c = 0..7:  Nx Ny Nxy  Mx My Mxy  Tx Ty
 *   N = membrane forces; M = moments with the convention M = D * kappa, kappa = (bx,x  by,y  bx,y+by,x),
 *       bx = ry, by = -rx (Batoz); T = transverse shear = k G t * gamma (gamma = w,x + bx, w,y + by)
 *       from the assumed shear strain field. Reliable when the element side does not exceed the
 *       thickness (h <= t); for h >> t it underestimates by up to 40% (in the Kirchhoff limit the
 *       side constraint only returns D*bs,ss). In that range the shear force should be recovered
 *       from the derivatives of the nodally averaged moments. */
int dkmshell_vertex_stresses(int nn, const double *x, const double *y, const DKMShellMaterial *m,
                             const double *ul, double *S);

/* Stress resultants at any point (natural coordinates xi, eta; for the triangle xi = L2, eta = L3). */
int dkmshell_point_stresses(int nn, const double *x, const double *y, const DKMShellMaterial *m,
                            const double *ul, double xi, double eta, double *S);

/* Element area */
double dkmshell_area(int nn, const double *x, const double *y);

#ifdef __cplusplus
}
#endif
#endif // dkmshellkernel_h

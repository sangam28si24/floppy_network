/* Units
 *   Positions  rx, ry      : reduced, in [0,1]     (actual = pos * L)
 *   Velocities vel[]       : actual units, [L/T]
 *   Bar lengths restLength : actual units, [L]  (L = sqrt(N), read from
 *                             the first line of the .dat file)
 *   Volume V                : L^2
 *   K:vv                    : [L/T^2] (acceleration)
 */

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/file.h>
#include <errno.h>

/* Compile-time constants */
#define N                  1600
#define DOUBLE_N           ((double)N)
#define DIM                2
#define X_COMP             0
#define Y_COMP             1
#define I_INDEX            0
#define J_INDEX            1

/* PBC reference length (positions stored in [0,1]) */
#define LENGTH             1.0
#define HALF_LENGTH        0.5

/* CG solver convergence tolerance */
#define CGSOLVER_TOL       1e-10

/* Position correction triggers when max |dr_alpha| / ell_alpha exceeds this.
 * Checks every step, so the network is never allowed to drift past this
 * tolerance unnoticed. */
#define POS_CORR_TOL       1e-8

#define HUGE_INT           10000000
#define UNI                ((double)rand() / ((double)RAND_MAX + 1.0))
#define TWO_PI             6.28318530717958647

/* Global simulation state */
double L;            /* actual box side length = sqrt(N), read from .dat file */
double invL;         /* 1/L */
double V;            /* L*L */
double strain;       /* accumulated shear strain (Lees-Edwards) */
double leDisplacement;

int contacts;        /* number of bars, Nc */

double *rx, *ry;     /* node positions, reduced units in [0,1] */
double *auxVec;      /* scratch buffer of size N*DIM, used by enOnVector */

/* Bond-network arrays, allocated in readNetwork */
int    *lookupBond;  /* lookupBond[2*l + I/J_INDEX] = node indices of bar l */
double *nij;         /* nij[DIM*l + X/Y_COMP] = unit vector along bar l */
double *restLength;  /* restLength[l], actual units = L * reduced_distance */
double *physVel = NULL;   /* size N*DIM, actual units [L/T] */

/* Verbosity flag for the CG solver; keep at 0 during production runs. */
static int verboseCG = 0;

/* Section 1: Network I/O */
int readNetwork(char *inFileName)
{
    int i, l;
    double dummy, dx, dy, r;
    FILE *inFile;

    rx      = (double *)malloc(sizeof(double) * N);
    ry      = (double *)malloc(sizeof(double) * N);
    auxVec  = (double *)malloc(sizeof(double) * N * DIM);

    inFile = fopen(inFileName, "rb");
    if (!inFile) {
        printf("couldn't open %s\n", inFileName);
        return 0;
    }

    (void)fscanf(inFile, "%lf %lf %lf", &L, &strain, &dummy);
    invL = 1.0 / L;
    V    = L * L;

    leDisplacement = fmod(strain, 1.0);
    if      (leDisplacement < -0.5) leDisplacement += 1.0;
    else if (leDisplacement >  0.5) leDisplacement -= 1.0;

    for (i = 0; i < N; i++)
        (void)fscanf(inFile, "%lf %lf %lf", &(rx[i]), &(ry[i]), &dummy);

    (void)fscanf(inFile, "%d %lf %lf", &contacts, &dummy, &dummy);

    /* Sanity check: the box side should be consistent with N if L = sqrt(N) convention is used. */
    {
        double L_expected = sqrt((double)N);
        if (fabs(L - L_expected) / L_expected > 0.01) {
            printf("WARNING: L=%.4g in file but sqrt(N)=%.4g (N=%d). "
                   "If your network has a different node count than the "
                   "compile-time N=%d, recompile with the correct N.\n",
                   L, L_expected, N, N);
        }
    }
    printf("Network: N=%d  L=%.4g  contacts=%d  z=%.4g\n",
           N, L, contacts, 2.0 * contacts / DOUBLE_N);

    lookupBond = (int    *)malloc(sizeof(int)    * 2 * contacts);
    nij        = (double *)malloc(sizeof(double) * contacts * DIM);
    restLength = (double *)malloc(sizeof(double) * contacts);

    for (l = 0; l < contacts; l++) {
        (void)fscanf(inFile, "%d %d %lf",
                     &(lookupBond[2*l + I_INDEX]),
                     &(lookupBond[2*l + J_INDEX]),
                     &(restLength[l]));
        i = lookupBond[2*l + I_INDEX];
        int j = lookupBond[2*l + J_INDEX];
        dx = rx[j] - rx[i];
        dy = ry[j] - ry[i];
        /* Lees-Edwards PBC (positions in [0,1]) */
        if      (dy >=  HALF_LENGTH) { dy -= LENGTH; dx -= leDisplacement; }
        else if (dy < -HALF_LENGTH)  { dy += LENGTH; dx += leDisplacement; }
        if      (dx >=  HALF_LENGTH)   dx -= LENGTH;
        else if (dx < -HALF_LENGTH)    dx += LENGTH;
        r = sqrt(dx*dx + dy*dy); 
        nij[DIM*l + X_COMP] = dx / r;
        nij[DIM*l + Y_COMP] = dy / r;
    }
    fclose(inFile);
    return 1;
}

void saveNetwork(char *networkFileName)
{
    int i, l;
    FILE *f = fopen(networkFileName, "wb");
    fprintf(f, "%.15g\t%.15g\t-1\n", L, strain);
    for (i = 0; i < N; i++)
        fprintf(f, "%.15g\t%.15g\t-1\n", rx[i], ry[i]);
    fprintf(f, "%d\t-1\t-1\n", contacts);
    for (l = 0; l < contacts; l++)
        fprintf(f, "%d\t%d\t%.15g\n",
                lookupBond[2*l+I_INDEX], lookupBond[2*l+J_INDEX], restLength[l]);
    fclose(f);
}

/* Frees everything readNetwork() allocated. Call once, at exit. */
void freeNetwork(void)
{
    free(rx); free(ry); free(auxVec);
    free(lookupBond); free(restLength);
    free(nij);
}

/* Section 2: Linear Operators S, S^T, N = SS^T
 *
 *   S       : contact-space <- node-space    (Nc x 2N)
 *   S^T     : node-space   <- contact-space  (2N x Nc)
 *   N = SS^T: contact-space <- contact-space (Nc x Nc)
 * All three act on the current nij[]; call updateGeometry() first if
 * positions have changed. */

/* res[Nc] = S|vec[2N]>, i.e. res_alpha = n_alpha . (v_j - v_i) */
void esOnVector(double *vec, double *res)
{
    int l, i, j;
    for (l = 0; l < contacts; l++) {
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        res[l] =  nij[DIM*l+X_COMP] * (vec[DIM*j+X_COMP] - vec[DIM*i+X_COMP])
                + nij[DIM*l+Y_COMP] * (vec[DIM*j+Y_COMP] - vec[DIM*i+Y_COMP]);
    }
}

/* res[2N] = S^T|vec[Nc]>, i.e. res_i -= f_alpha n_alpha, res_j += f_alpha n_alpha.
 * This is the adjoint of esOnVector: it scatters contact-space forces back onto the nodes. */
void esTransposeOnVector(double *vec, double *res)
{
    int l, i, j;
    for (i = 0; i < N*DIM; i++) res[i] = 0.0;
    for (l = 0; l < contacts; l++) {
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        res[DIM*i + X_COMP] -= vec[l] * nij[DIM*l + X_COMP];
        res[DIM*j + X_COMP] += vec[l] * nij[DIM*l + X_COMP];
        res[DIM*i + Y_COMP] -= vec[l] * nij[DIM*l + Y_COMP];
        res[DIM*j + Y_COMP] += vec[l] * nij[DIM*l + Y_COMP];
    }
}

/* res[Nc] = N|vec[Nc]> = SS^T|vec>, via the auxVec scratch buffer*/
void enOnVector(double *vec, double *res)
{
    esTransposeOnVector(vec, auxVec);
    esOnVector(auxVec, res);
}

/* Contact-space dot product, used only inside cgSolver. */
double cdot(double *v1, double *v2)
{
    int l;
    double s = 0.0;
    for (l = 0; l < contacts; l++) s += v1[l] * v2[l];
    return s;
}

/* Conjugate-gradient solve of N|x> = |b>, no preconditioner.
 * x must already be initialized (pass a zeroed buffer for a zero start).
 * Convergence is judged on the max-norm residual relative to ||b||.
 * If verboseCG is set, progress prints to stderr on one overwritten line.
 * Returns 1 on convergence, 0 if it gives up after HUGE_INT iterations. */
int cgSolver(double *x, double *b)
{
    double *r, *w, *z;
    double a, B, maxNorm, scale;
    int i, k;

    if (contacts == 0) return 1;   /* nothing to solve */

    scale = sqrt(cdot(b, b) / (double)contacts);
    if (scale < 1e-30) return 1;   /* b is effectively zero, trivial solution */

    r = (double *)malloc(contacts * sizeof(double));
    w = (double *)malloc(contacts * sizeof(double));
    z = (double *)malloc(contacts * sizeof(double));

    enOnVector(x, r);
    for (i = 0; i < contacts; i++) { r[i] = b[i] - r[i]; w[i] = -r[i]; }

    enOnVector(w, z);
    a = cdot(r, w) / cdot(w, z);
    for (i = 0; i < contacts; i++) x[i] += a * w[i];

    for (k = 0; k < HUGE_INT; k++) {
        maxNorm = 0.0;
        for (i = 0; i < contacts; i++) {
            r[i] -= a * z[i];
            if (r[i]*r[i] > maxNorm) maxNorm = r[i]*r[i];
        }
        maxNorm = sqrt(maxNorm);
        if (verboseCG)
            fprintf(stderr, "\r  CG iter %6d  res/scale = %.4e  ", k, maxNorm/scale);
        if (maxNorm / scale < CGSOLVER_TOL) break;

        B = cdot(r, z) / cdot(w, z);
        for (i = 0; i < contacts; i++) w[i] = -r[i] + B * w[i];
        enOnVector(w, z);
        a = cdot(r, w) / cdot(w, z);
        for (i = 0; i < contacts; i++) x[i] += a * w[i];
    }
    if (verboseCG) fprintf(stderr, "\n");
    free(r); free(w); free(z);
    if (k == HUGE_INT) { printf("WARNING: CG did not converge!\n"); return 0; }
    return 1;
}

/* Convenience wrapper: allocate a zeroed contact-space buffer, solve
 * N f = b, return f. Caller owns the returned pointer and must free it. */
static double *solveN(double *b)
{
    double *f = (double *)calloc(contacts, sizeof(double));
    cgSolver(f, b);
    return f;
}

/* Section 3: Geometry Update */

/*Recomputes nij[] from the current rx[], ry[]. Must be called after any position change (drift step, position correction, affine shear step)
 -- every linear operator above implicitly depends on nij, so
 * stale unit vectors silently produce wrong dynamics.
 * Uses the same Lees-Edwards PBC convention as readNetwork()*/
void updateGeometry(void)
{
    int l, i, j;
    double dx, dy, r;
    for (l = 0; l < contacts; l++) {
        i  = lookupBond[2*l + I_INDEX];
        j  = lookupBond[2*l + J_INDEX];
        dx = rx[j] - rx[i];
        dy = ry[j] - ry[i];
        if      (dy >=  HALF_LENGTH) { dy -= LENGTH; dx -= leDisplacement; }
        else if (dy < -HALF_LENGTH)  { dy += LENGTH; dx += leDisplacement; }
        if      (dx >=  HALF_LENGTH)   dx -= LENGTH;
        else if (dx < -HALF_LENGTH)    dx += LENGTH;
        r = sqrt(dx*dx + dy*dy);
        nij[DIM*l + X_COMP] = dx / r;
        nij[DIM*l + Y_COMP] = dy / r;
    }
}

/* Computes actual bar lengths: lengths[l] = L * sqrt(dx_reduced^2 + dy_reduced^2).
 * Used by maxBarLengthError() and applyPositionCorrection(). */
static void computeBarLengths(double *lengths)
{
    int l, i, j;
    double dx, dy;
    for (l = 0; l < contacts; l++) {
        i  = lookupBond[2*l + I_INDEX];
        j  = lookupBond[2*l + J_INDEX];
        dx = rx[j] - rx[i];
        dy = ry[j] - ry[i];
        if      (dy >=  HALF_LENGTH) { dy -= LENGTH; dx -= leDisplacement; }
        else if (dy < -HALF_LENGTH)  { dy += LENGTH; dx += leDisplacement; }
        if      (dx >=  HALF_LENGTH)   dx -= LENGTH;
        else if (dx < -HALF_LENGTH)    dx += LENGTH;
        lengths[l] = L * sqrt(dx*dx + dy*dy);   /* actual units */
    }
}

/* Section 4: Dynamics Building Blocks */

/* Kinetic energy, KE = (1/2) * sum_i |v_i|^2 (unit mass, actual velocity units). */
double kineticEnergy(double *vel)
{
    int k;
    double ke = 0.0;
    for (k = 0; k < N*DIM; k++) ke += vel[k] * vel[k];
    return 0.5 * ke;
}

/* Total system momentum, |P> = sum_i m_i |v_i> = sum_i |v_i> (unit mass).
 * Should sit at zero to floating-point precision for the whole run: the
 * initial condition has zero COM velocity by construction (see
 * initializePhysicalVelocity), and every internal force in this model
 * comes from S^T|f>, which places equal-and-opposite contributions on
 * each bar's two endpoints (Newton's third law is built into S^T's
 * definition) -- so RATTLE never has a mechanism to inject net momentum.
 * This diagnostic exists to verify that expectation empirically rather
 * than assume it. */
void computeMomentum(double *vel, double *px, double *py)
{
    int i;
    double sx = 0.0, sy = 0.0;
    for (i = 0; i < N; i++) {
        sx += vel[DIM*i + X_COMP];
        sy += vel[DIM*i + Y_COMP];
    }
    *px = sx;
    *py = sy;
}

/* Curvature contraction (K:vv)_alpha for every bar:
 *
 *   (K:vv)_alpha = (1/r_alpha) * [ |dv|^2 - (dv . n_alpha)^2 ]
 *                = (1/r_alpha) * |dv_perp|^2
 *
 * where dv = vel[j] - vel[i] and r_alpha is approximated by
 * restLength[l] (cheaper than recomputing the actual current length,
 * and accurate since the network stays close to its rest lengths).
 * Physically: this is how much a bond's rigid-body rotation "steers"
 * it off a straight line as the nodes move, so a satisfied bar needs a
 * centripetal tension to stay at fixed length. Result is in units of
 * acceleration, [L/T^2]. See README.md for the full derivation. */
void computeKvv(double *vel, double *kvv)
{
    int l, i, j;
    double dux, duy, un, perp_sq;
    for (l = 0; l < contacts; l++) {
        i   = lookupBond[2*l + I_INDEX];
        j   = lookupBond[2*l + J_INDEX];
        dux = vel[DIM*j + X_COMP] - vel[DIM*i + X_COMP];
        duy = vel[DIM*j + Y_COMP] - vel[DIM*i + Y_COMP];
        un  = nij[DIM*l+X_COMP]*dux + nij[DIM*l+Y_COMP]*duy; /* component of dv along the bond */
        perp_sq = dux*dux + duy*duy - un*un;
        kvv[l] = perp_sq / restLength[l];
    }
}

/* Node-space acceleration |a> = -S^T N^{-1} |K:vv>: the centripetal
 * force field needed to keep every bar at its rest length given the
 * current velocities. Steps: solve N f = K:vv in contact space, then
 * scatter f back to nodes with -S^T. This is the single most expensive
 * call in the program (one CG solve per call), and it's called at
 * least once per RATTLE step. */
void computeAcceleration(double *vel, double *acc)
{
    double *kvv = (double *)malloc(contacts * sizeof(double));
    int k;

    computeKvv(vel, kvv);
    double *f = solveN(kvv);
    esTransposeOnVector(f, acc);
    for (k = 0; k < N*DIM; k++) acc[k] = -acc[k];

    free(kvv);
    free(f);
}

/* Projects vel onto the constraint manifold ker(S):
 *   |v> <- |v> - S^T N^{-1} S|v>
 * Removes whatever component of vel violates S|v> = 0 (i.e. would
 * stretch or compress a bar). After this call, esOnVector(vel, ...)
 * should return numerically zero for every bar. */
void projectVelocities(double *vel)
{
    double *sv   = (double *)malloc(contacts * sizeof(double));
    double *corr = (double *)malloc(N*DIM    * sizeof(double));
    int k;

    esOnVector(vel, sv);          /* sv = S|vel> */
    double *f = solveN(sv);       /* f  = N^{-1} S|vel> */
    esTransposeOnVector(f, corr); /* corr = S^T f = S^T N^{-1} S|vel> */
    for (k = 0; k < N*DIM; k++) vel[k] -= corr[k];

    free(sv); free(f); free(corr);
}

/* Max relative bar-length error, max_alpha |r_alpha - ell_alpha| / ell_alpha.
 * Diagnostic for how far the simulation has drifted off the constraint manifold. */
double maxBarLengthError(void)
{
    double *rCur = (double *)malloc(contacts * sizeof(double));
    double maxErr = 0.0, err;
    int l;
    computeBarLengths(rCur);
    for (l = 0; l < contacts; l++) {
        err = fabs(rCur[l] - restLength[l]) / restLength[l];
        if (err > maxErr) maxErr = err;
    }
    free(rCur);
    return maxErr;
}

/* Max velocity constraint violation, max_alpha |(S|v>)_alpha|. 
 * Should sit at machine epsilon (~1e-15) right after every RATTLE velocity projection; a larger value usually means something upstream forgot to project. */
double maxVelocityConstraint(double *vel)
{
    double *sv = (double *)malloc(contacts * sizeof(double));
    double maxV = 0.0, av;
    int l;
    esOnVector(vel, sv);
    for (l = 0; l < contacts; l++) {
        av = fabs(sv[l]);
        if (av > maxV) maxV = av;
    }
    free(sv);
    return maxV;
}

/* One-pass position correction, applied whenever the bar-length drift
 * that RATTLE leaves behind (O(dt^2) per step) crosses POS_CORR_TOL.
 * The kinetic energy is recorded first, then the current bar-length
 * errors dr_alpha = r_alpha(x) - ell_alpha are used to solve N f = dr
 * for the minimum-norm correction dX = -S^T f. Positions are shifted by
 * dX, the geometry is rebuilt, velocities are re-projected onto the
 * corrected manifold, and finally rescaled to restore KE to its value
 * from before the correction. Returns the max relative bar-length error
 * measured before the correction, so callers can log how bad the drift
 * had gotten. */
double applyPositionCorrection(double *vel)
{
    double *rCur = (double *)malloc(contacts * sizeof(double));
    double *dr   = (double *)malloc(contacts * sizeof(double));
    double *dX   = (double *)malloc(N*DIM    * sizeof(double));
    int l, k, i;

    /* Step 1: record kinetic energy */
    double E0 = kineticEnergy(vel);

    /* Step 2: bar-length errors dr_alpha (actual units) */
    computeBarLengths(rCur);
    double maxErr = 0.0;
    for (l = 0; l < contacts; l++) {
        dr[l] = rCur[l] - restLength[l];
        double relErr = fabs(dr[l]) / restLength[l];
        if (relErr > maxErr) maxErr = relErr;
    }
    free(rCur);

    /* Step 3: solve N f = dr, compute dX = -S^T f */
    double *f = solveN(dr);
    esTransposeOnVector(f, dX);
    free(f); free(dr);

    /* Step 4: update positions (reduced) */
    for (i = 0; i < N; i++) {
        double newrx = rx[i] - dX[DIM*i + X_COMP] * invL;
        double newry = ry[i] - dX[DIM*i + Y_COMP] * invL;
        /* wrap back into [0,1] */
        if      (newrx >= LENGTH) newrx -= LENGTH;
        else if (newrx <  0.0)   newrx += LENGTH;
        if      (newry >= LENGTH) newry -= LENGTH;
        else if (newry <  0.0)   newry += LENGTH;
        rx[i] = newrx;
        ry[i] = newry;
    }
    free(dX);

    /* Step 5: rebuild geometry */
    updateGeometry();

    /* Step 6: re-project velocities */
    projectVelocities(vel);

    /* Step 7: rescale to restore energy (projection can only reduce KE, so lambda >= 1) */
    double KEnew = kineticEnergy(vel);
    if (KEnew > 1e-30) {
        double lambda = sqrt(E0 / KEnew);
        for (k = 0; k < N*DIM; k++) vel[k] *= lambda;
    }

    return maxErr;
}

/* One full RATTLE step (velocity-Verlet, symplectic). The half-kick
 * v^{n+1/2} = v^n + (dt/2) a^n is followed by the drift x^{n+1} = x^n +
 * dt * v^{n+1/2} (in reduced units) and a rebuild of the geometry at
 * x^{n+1}. The new acceleration a^{n+1} is evaluated there using the
 * mid-step velocity v^{n+1/2}, then the second half-kick brings the
 * velocity up to v^{n+1} = v^{n+1/2} + (dt/2) a^{n+1}, and a final
 * projection v^{n+1} <- v^{n+1} - S^T N^{-1} S v^{n+1} removes the
 * O(dt^2) drift the drift step introduced off the constraint manifold.
 * On entry: vel = v^n, acc = a^n. On return: vel = v^{n+1}, acc = a^{n+1}. */
void rattleStep(double *vel, double *acc, double dt)
{
    int i, k;
    double half_dt = 0.5 * dt;

    /* 1. Half-kick */
    for (k = 0; k < N*DIM; k++) vel[k] += half_dt * acc[k];

    /* 2. Drift: update reduced positions using actual velocities / L */
    for (i = 0; i < N; i++) {
        double newrx = rx[i] + dt * vel[DIM*i + X_COMP] * invL;
        double newry = ry[i] + dt * vel[DIM*i + Y_COMP] * invL;
        if      (newrx >= LENGTH) newrx -= LENGTH;
        else if (newrx <  0.0)   newrx += LENGTH;
        if      (newry >= LENGTH) newry -= LENGTH;
        else if (newry <  0.0)   newry += LENGTH;
        rx[i] = newrx;
        ry[i] = newry;
    }

    /* 3. Rebuild geometry at x^{n+1} */
    updateGeometry();

    /* 4. Centripetal acceleration using half-step velocities at the new geometry */
    computeAcceleration(vel, acc);   /* vel is still v^{n+1/2} here */

    /* 5. Second half-kick */
    for (k = 0; k < N*DIM; k++) vel[k] += half_dt * acc[k];

    /* 6. Remove the O(dt^2) velocity drift off the constraint manifold */
    projectVelocities(vel);
}

/* Instantaneous shear stress,
 *   sigma_xy = (1/V) * sum_alpha (N^{-1} K:vv)_alpha * n_{alpha,x} * n_{alpha,y} * ell_alpha,
 * the virial contribution from the rigid-bar tensions. Note that
 * (N^{-1} K:vv)_alpha is the negative of the bar tension in
 * computeCurvatureForce()'s convention (tension = -N^{-1} K:vv); the two
 * sign conventions cancel, so this returns the same physical stress as
 * computeCurvatureStress(). */
double instantaneousStress(double *vel)
{
    double *kvv = (double *)malloc(contacts * sizeof(double));
    computeKvv(vel, kvv);
    double *f = solveN(kvv);
    free(kvv);

    double sigma = 0.0;
    int l;
    for (l = 0; l < contacts; l++)
        sigma += f[l] * nij[DIM*l+X_COMP] * nij[DIM*l+Y_COMP] * restLength[l];
    free(f);
    return sigma / V;
}

/* Random floppy-space velocity initialization.
 * Draws a Gaussian node velocity field (Box-Muller, sigma = sqrt(T)),
 * removes the center-of-mass drift, projects onto the floppy (constraint)
 * manifold via
 *   |Rdot> = P|V_random> = (I - S^T N^{-1} S) |V_random>,   N = SS^T
 * and rescales so that KE = (DOF/2) * T, with
 *   DOF = 2N - contacts - DIM
 * The -DIM accounts for the d=2 rigid-body translation (COM) modes removed
 * below.
 *
 * COM removal is necessary and must be explicit: COM translation lives in
 * ker(S) (it changes no bar length), so the floppy projector P *retains* it
 * rather than removing it. Without explicit subtraction, the whole network
 * drifts uniformly, contributing KE without any physical thermal content,
 * and biasing the temperature rescale.
 *
 * After the call, vel satisfies S|vel> = 0 and has zero COM velocity. */
void initializeVelocities(double *vel, double temperature)
{
    int i, k;
    double sigma = sqrt(temperature);

    /* Box-Muller Gaussian samples: raw |V_random> */
    for (i = 0; i < N; i++) {
        double u1 = UNI + 1e-15;
        double u2 = UNI;
        double r  = sigma * sqrt(-2.0 * log(u1));
        vel[DIM*i + X_COMP] = r * cos(TWO_PI * u2);
        vel[DIM*i + Y_COMP] = r * sin(TWO_PI * u2);
    }

    /* Remove center-of-mass velocity before projection.
     * COM translation is in ker(S) so P would leave it intact;
     * subtract it explicitly so the network has zero net drift. */
    double com_x = 0.0, com_y = 0.0;
    for (i = 0; i < N; i++) {
        com_x += vel[DIM*i + X_COMP];
        com_y += vel[DIM*i + Y_COMP];
    }
    com_x /= (double)N;
    com_y /= (double)N;
    for (i = 0; i < N; i++) {
        vel[DIM*i + X_COMP] -= com_x;
        vel[DIM*i + Y_COMP] -= com_y;
    }

    /* Project onto the constraint manifold: |Rdot> = (I - S^T N^{-1} S) |V_random> */
    projectVelocities(vel);

    /* Rescale to target temperature. In DOF = 2N - contacts - DIM, 2N is
     * the total node-space dimensions, contacts is subtracted because
     * each rigid-bar constraint removes one mode, and DIM is subtracted
     * for the d center-of-mass translation modes removed above; then
     * <KE> = (DOF/2) * T by equipartition. */
    int DOF = 2*N - contacts - DIM;
    if (DOF <= 0) DOF = 1;            /* safety: rigid or over-constrained network */
    double KE_target = 0.5 * DOF * temperature;
    double KE_cur    = kineticEnergy(vel);
    if (KE_cur > 1e-30) {
        double scale = sqrt(KE_target / KE_cur);
        for (k = 0; k < N*DIM; k++) vel[k] *= scale;
    }
}

/* Initializes the single persistent physical velocity, exactly once. The
 * first call allocates physVel and draws it via initializeVelocities()
 * above (Box-Muller, projection, DOF/2*T rescale); every later call is a
 * deliberate no-op, so runNVE() and runGmod() can each just ask for
 * physVel to be ready without risking an accidental re-draw. */
void initializePhysicalVelocity(double temperature)
{
    if (physVel != NULL) return;    /* already initialized, never redraw */
    physVel = (double *)malloc(N*DIM * sizeof(double));
    initializeVelocities(physVel, temperature);
}

/* Section 5: Nonlinear (Curvature-Based) Shear Modulus
 *
 * A floppy network (z < zc) has zero linear elastic response: to first
 * order in an imposed affine shear strain gamma, the bond-length change
 * can always be cancelled by a compatible non-affine velocity field,
 * which is exactly the floppy projector already used in
 * projectVelocities(). The obstruction shows up at second order. Moving
 * nodes along the current physical velocity V for a finite step bends
 * each bond's direction n_ij, so a rigid bond of fixed length sweeps out
 * a curved path, i.e. (K:VV)_alpha != 0 even though S|V> = 0 exactly.
 * Since bars can't stretch, an internal bar tension |f> is required to
 * exactly cancel that curvature, a centripetal Lagrange multiplier,
 *
 *     N f = -(K:VV),   N = SS^T,   f = -N^{-1}|K:VV>,
 *
 * the same computeKvv() / N^{-1} pairing computeAcceleration() uses for
 * the RATTLE dynamics, just evaluated for its own sake here. The
 * resulting bond tension is contracted against the affine strain-rate
 * vector s_alpha = n_x n_y ell_alpha to get the curvature stress,
 *
 *     Gamma = Gamma(R, V) = -<s | f> / V,
 *
 * evaluated directly from the current physical velocity, with no
 * auxiliary probe velocity drawn for this. The nonlinear shear modulus
 * is the genuine two-state finite difference: measure Gamma1 =
 * Gamma(R,V), apply a finite affine shear step, restore exact bond
 * lengths with the iterated Lerner correction
 * (relaxToConstraintManifold), re-project the same physical velocity V
 * onto the updated rigidity matrix with no new velocity drawn, measure
 * Gamma2 = Gamma(R',V), and difference,
 *
 *     G = (Gamma2 - Gamma1) / dgamma.
 *
 * Because V is the single persistent physical velocity carried through
 * the whole run, computeShearModulusCurvature() below permanently
 * mutates both the geometry (sheared and relaxed) and V (re-projected):
 * it is a real physical step, not a side-effect-free probe. Callers
 * comparing several dgamma values from the same base state, like
 * runGmod's sweep, need to snapshot and restore themselves between
 * trials.
 */

/* s_alpha = dr_alpha/dgamma at the affine level = n_x n_y ell_alpha.
 * This is the affine-shear rate vector fed into the curvature stress. */
void computeAffineStrainVector(double *s)
{
    int l;
    for (l = 0; l < contacts; l++)
        s[l] = nij[DIM*l+X_COMP] * nij[DIM*l+Y_COMP] * restLength[l];
}

/* Curvature-based bond tension |f> = -N^{-1}|K:VV>. Vfield must already
 * be floppy-projected (S|Vfield> = 0), which is always true for the
 * persistent physical velocity: it's projected at initialization and
 * re-projected after every RATTLE step or geometric change, never used
 * raw. */
void computeCurvatureForce(double *Vfield, double *f)
{
    double *kvv = (double *)malloc(contacts * sizeof(double));
    double *tmp;
    int l;

    computeKvv(Vfield, kvv);      /* (K:VV)_alpha */
    tmp = solveN(kvv);            /* tmp = N^{-1}|K:VV>, N = SS^T */
    for (l = 0; l < contacts; l++) f[l] = -tmp[l];   /* f = -N^{-1}|K:VV> */

    free(kvv); free(tmp);
}

/* Curvature stress Gamma(R,V), evaluated directly from the current
 * velocity `vel`, which is already on the constraint manifold, so no
 * projection or random draw happens here. The bond tension
 * |f> = -N^{-1}|K:VV(vel)> is contracted against the affine strain-rate
 * vector s_alpha = n_x n_y ell_alpha to give Gamma = -<s|f> / V. */
double computeCurvatureStress(double *vel)
{
    double *s = (double *)malloc(contacts * sizeof(double));
    double *f = (double *)malloc(contacts * sizeof(double));
    double Gamma;

    computeCurvatureForce(vel, f);
    computeAffineStrainVector(s);
    Gamma = -cdot(s, f) / V;

    free(s); free(f);
    return Gamma;
}

/* Pressure p(R,V) -- the isotropic twin of Gamma(R,V) above, built from
 * the exact same bond tension |f> = -N^{-1}|K:VV(vel)>, just contracted
 * against a different kinematic vector.
 *
 * Gamma uses s_alpha = dr_alpha/dgamma under simple shear (n_x n_y ell_alpha).
 * Pressure needs dr_alpha/depsilon under isotropic compression instead:
 * with x_i -> (1-epsilon) x_i, bond alpha=(i,j) scales the same way, so
 * r_alpha(epsilon) = (1-epsilon) r_alpha(0), giving dr_alpha/depsilon = -r_alpha.
 * Substituting into the same conversion rule used for Gamma
 * (Sigma_eta = -<ds/deta|f>/V) gives Sigma_epsilon = -<-r|f>/V = +<r|f>/V.
 * Dividing by dimension DIM for the trace-averaged convention gives:
 *
 *     p = <r|f> / (V * DIM) = (1/(V*DIM)) sum_alpha r_alpha f_alpha
 *
 * This matches [LDW]'s own pressure formula (their Eq. 24,
 * p = <f|r>/(Omega*d), no leading minus) even though Gamma needed one:
 * the difference isn't a separate convention, it's that s_epsilon = -r
 * already carries an intrinsic minus (contraction kinematics) that
 * cancels the same -<.|f>/V rule's minus, while s_alpha = +n_x n_y ell_alpha
 * (shear kinematics) does not. Same rule, two different kinematic signs.
 *
 * r_alpha (current bond length) comes from computeBarLengths(), already
 * used by maxBarLengthError()/rmsBarLengthError(); f comes from
 * computeCurvatureForce(), already used by computeCurvatureStress() --
 * nothing new is computed here, just a different dot product. */
double computePressure(double *vel)
{
    double *r = (double *)malloc(contacts * sizeof(double));
    double *f = (double *)malloc(contacts * sizeof(double));
    double p;

    computeBarLengths(r);
    computeCurvatureForce(vel, f);
    p = cdot(r, f) / (V * (double)DIM);

    free(r); free(f);
    return p;
}


/* Applies a finite affine simple-shear increment dgamma to every node
 * position, x_i <- x_i + dgamma * y_i, and updates the strain
 * bookkeeping so later updateGeometry() calls compute nij consistently.
 * Positions are stored in reduced units [0,1), so the affine shift in
 * reduced units is drx_i = dgamma * ry[i]: the actual shift is
 * dgamma * Y_i, and reduced = actual / L. */
void applyAffineShearStep(double dgamma)
{
    int i;
    for (i = 0; i < N; i++) {
        double newrx = rx[i] + dgamma * ry[i];
        if      (newrx >= LENGTH) newrx -= LENGTH;
        else if (newrx <  0.0)    newrx += LENGTH;
        rx[i] = newrx;
    }
    strain += dgamma;
    leDisplacement = fmod(strain, 1.0);
    if      (leDisplacement < -0.5) leDisplacement += 1.0;
    else if (leDisplacement >  0.5) leDisplacement -= 1.0;
    updateGeometry();
}

#define RELAX_MAX_ITERS   200
#define RELAX_TOL         1e-12

/* Iteratively enforces the exact rigid-bar constraint r_alpha(x) =
 * ell_alpha for every bar, via repeated minimum-norm corrections
 *     dX = -S^T N^{-1} |dr>
 * -- the same displacement field applyPositionCorrection() uses once
 * per call for drift cleanup during NVE, but here iterated to full
 * Newton convergence instead of just one pass. Purely geometric: no
 * velocities are touched.
 *
 * Returns the number of iterations used to reach RELAX_TOL, or -1 if
 * it didn't converge within RELAX_MAX_ITERS. A -1 is itself a
 * meaningful result: it means the network ran out of floppy modes and
 * couldn't accommodate the imposed strain non-affinely, i.e. it went
 * isostatic or self-stressed at this strain level. */
int relaxToConstraintManifold(void)
{
    double *rCur = (double *)malloc(contacts * sizeof(double));
    double *dr   = (double *)malloc(contacts * sizeof(double));
    double *dX   = (double *)malloc(N*DIM    * sizeof(double));
    int it, l, i;
    int result = -1;

    for (it = 0; it < RELAX_MAX_ITERS; it++) {
        double maxErr = 0.0;
        double *f;

        computeBarLengths(rCur);
        for (l = 0; l < contacts; l++) {
            dr[l] = rCur[l] - restLength[l];
            double relErr = fabs(dr[l]) / restLength[l];
            if (relErr > maxErr) maxErr = relErr;
        }
        if (maxErr < RELAX_TOL) { result = it; break; }

        f = solveN(dr);
        esTransposeOnVector(f, dX);
        free(f);

        for (i = 0; i < N; i++) {
            double newrx = rx[i] - dX[DIM*i + X_COMP] * invL;
            double newry = ry[i] - dX[DIM*i + Y_COMP] * invL;
            if      (newrx >= LENGTH) newrx -= LENGTH;
            else if (newrx <  0.0)   newrx += LENGTH;
            if      (newry >= LENGTH) newry -= LENGTH;
            else if (newry <  0.0)   newry += LENGTH;
            rx[i] = newrx; ry[i] = newry;
        }
        updateGeometry();
    }

    free(rCur); free(dr); free(dX);
    return result;
}

/* The genuine nonlinear (curvature-based) shear modulus:
 *
 *     Gamma1 = Gamma(R, V)                      [current physical V]
 *     x <- x + dgamma * y                       [affine shear]
 *     Lerner correction, iterated to tolerance   [relaxToConstraintManifold]
 *     rebuild S                                  [automatic, via updateGeometry]
 *     V <- (I - S^T N^-1 S) V                    [same V, re-projected]
 *     Gamma2 = Gamma(R', V)
 *     G = (Gamma2 - Gamma1) / dgamma
 *
 * `vel` is the single persistent physical velocity for the whole run
 * (see initializePhysicalVelocity / physVel), passed by pointer and
 * genuinely mutated in place by the re-projection step -- no new
 * random velocity is drawn here. The geometry (rx, ry, strain) is
 * likewise left in the sheared, Lerner-relaxed state: this is a real
 * physical step, not a side-effect-free probe. Callers wanting several
 * dgamma trials from the same base state (e.g. a convergence sweep)
 * need to snapshot/restore (rx, ry, strain, leDisplacement, vel)
 * themselves around each call.
 *
 * *relax_iters_out (if non-NULL) receives the number of Newton
 * iterations relaxToConstraintManifold() needed, or -1 on failure.
 * *residual_out (if non-NULL) receives the max relative bar-length
 * error at the strained, relaxed state -- the physically meaningful
 * convergence diagnostic for this dgamma. */
double computeShearModulusCurvature(double *vel, double dgamma,
                                     int *relax_iters_out,
                                     double *residual_out)
{
    double Gamma1 = computeCurvatureStress(vel);   /* Gamma(R,V) */

    applyAffineShearStep(dgamma);                  /* x <- x + dgamma*y */
    int iters = relaxToConstraintManifold();        /* Lerner correction, iterated;
                                                         also rebuilds S via updateGeometry */

    projectVelocities(vel);                         /* V <- (I - S^T N^-1 S) V, same V */

    double Gamma2 = computeCurvatureStress(vel);    /* Gamma(R',V) */

    if (residual_out)    *residual_out    = maxBarLengthError();
    if (relax_iters_out) *relax_iters_out = iters;

    return (Gamma2 - Gamma1) / dgamma;
}

/* Section 6: NVE Run */

/*
 * Runs RATTLE NVE dynamics and writes diagnostics to file.
 *
 * Operates on the single persistent physical velocity `physVel`:
 * initializePhysicalVelocity() draws it once (no-op if an earlier call
 * in this run already did), then every RATTLE step evolves that same
 * buffer in place -- vel is never reallocated or redrawn here.
 *
 * temperature sets the initial temperature and is only used if physVel
 * hasn't already been initialized, dt is the timestep in actual units
 * (typically 0.001-0.005), N_steps is the total number of RATTLE steps,
 * out_every controls how often diagnostics are written, and outFileName
 * is the path for the output data file. Each output row records step,
 * KE, max_bar_err, max_vel_constraint, sigma_xy, pressure, px, py,
 * tab-separated. pressure (computePressure) is the isotropic twin of
 * sigma_xy -- same bond tension, dotted against bond length instead of
 * the shear kinematic vector -- so pressure-vs-stress is just two
 * columns of the same file. px, py (computeMomentum) are the total
 * system momentum -- should stay at ~0 for the whole run; see
 * computeMomentum()'s comment for why.
 */
void runNVE(double temperature, double dt, int N_steps,
            int out_every, char *outFileName)
{
    int step;
    double *acc = (double *)malloc(N*DIM * sizeof(double));
    FILE   *out;

    printf("\n=== NVE RATTLE dynamics ===\n");
    printf("  T=%.4g  dt=%.4g  steps=%d  write_every=%d\n",
           temperature, dt, N_steps, out_every);
    printf("  DOF = 2N - Nc - d = %d\n", 2*N - contacts - DIM);
    printf("  Output: %s\n\n", outFileName);

    /* Initialize the single persistent physical velocity (no-op if
     * already initialized by an earlier call in this run) */
    initializePhysicalVelocity(temperature);
    double *vel = physVel;
    computeAcceleration(vel, acc);

    /* Open output */
    out = fopen(outFileName, "w");
    fprintf(out, "# rigid-bar NVE RATTLE  N=%d  contacts=%d  z=%.4g\n",
            N, contacts, 2.0*contacts/DOUBLE_N);
    fprintf(out, "# T=%.6g  dt=%.6g\n", temperature, dt);
    fprintf(out, "# col: step  KE  max_bar_err  max_vel_constraint  sigma_xy  pressure  px  py\n");

    /* Main loop */
    for (step = 0; step <= N_steps; step++) {

        /* diagnostics & output */
        if (step % out_every == 0) {
            double KE       = kineticEnergy(vel);
            double barErr   = maxBarLengthError();
            double velCon   = maxVelocityConstraint(vel);
            double sigma    = instantaneousStress(vel);
            double pressure = computePressure(vel);
            double px, py;
            computeMomentum(vel, &px, &py);
            fprintf(out, "%d\t%.10g\t%.6e\t%.6e\t%.10g\t%.10g\t%.6e\t%.6e\n",
                    step, KE, barErr, velCon, sigma, pressure, px, py);
            fflush(out);
            printf("  step %7d  KE=%.6g  barErr=%.2e  velCon=%.2e  sigma=%.4g  "
                   "pressure=%.4g  |P|=%.2e\n",
                   step, KE, barErr, velCon, sigma, pressure, sqrt(px*px + py*py));
        }

        if (step == N_steps) break;

        /* Position correction, checked every step. maxBarLengthError() is
         * just a loop over bonds (no linear solve), so checking it here
         * costs nothing next to a RATTLE step; applyPositionCorrection()
         * itself only runs when the drift has actually crossed the
         * tolerance. */
        if (maxBarLengthError() > POS_CORR_TOL) {
            applyPositionCorrection(vel);
            computeAcceleration(vel, acc);
        }

        /* RATTLE step: vel (== physVel) evolves in place, V^n -> V^n+1 */
        rattleStep(vel, acc, dt);
    }

    fclose(out);
    printf("\nNVE run complete. Output in %s\n", outFileName);

    /* physVel is intentionally NOT freed here -- it's the single
     * persistent physical velocity for the whole run and may still be
     * needed afterward (e.g. a gmod measurement right after). It gets
     * freed exactly once, in main(), at the very end of the program. */
    free(acc);
}

/* Section 7: Nonlinear G(z) Driver */

/* Advances the single persistent physical velocity by N_steps of RATTLE
 * dynamics, applying the same position correction runNVE() uses, but
 * with no diagnostic file output. This is the thermalization step the
 * curvature-modulus measurement needs before it can measure anything
 * meaningful: G(z) should reflect a dynamically evolved state (R,V), not
 * a velocity field that was just drawn and never allowed to interact
 * with the geometry. It reuses rattleStep() and applyPositionCorrection()
 * exactly as runNVE() does, evolving `vel` (== physVel) and the geometry
 * in place for N_steps before returning control to the caller for
 * measurement. */
void advanceDynamicsToMeasurementState(double *vel, double dt, int N_steps)
{
    double *acc = (double *)malloc(N*DIM * sizeof(double));
    int step;

    computeAcceleration(vel, acc);
    for (step = 0; step < N_steps; step++) {
        if (maxBarLengthError() > POS_CORR_TOL) {
            applyPositionCorrection(vel);
            computeAcceleration(vel, acc);
        }
        rattleStep(vel, acc, dt);
    }
    free(acc);
}

/* Same dynamics loop as advanceDynamicsToMeasurementState() above --
 * identical physics, same per-step drift correction -- but additionally
 * samples Gamma = instantaneousStress(vel) periodically during the
 * SECOND HALF of the window (the first half is discarded as a
 * transient, letting velocity settle into the new geometry before
 * anything is averaged) and reports the sample mean and standard error.
 *
 * This exists to test directly whether a single endpoint snapshot of
 * Gamma is simply too noisy, on its own, to reveal a strain-dependent
 * trend -- since Gamma is purely kinetic (built from K:vv), a lone
 * instantaneous reading can fluctuate with the sign and size of
 * whatever the velocity happens to be doing at that moment, the same
 * way one gas-molecule snapshot doesn't show you the pressure of a gas.
 * Averaging over many samples along the trajectory is the direct way
 * to check whether a real trend is hiding underneath that noise.
 *
 * sample_every controls how many RATTLE steps separate consecutive
 * samples (use >1 to reduce serial correlation between samples and cut
 * the cost of the instantaneousStress() solves; 1 samples every step). */
void advanceDynamicsAndAverageStress(double *vel, double dt, int N_steps,
                                      int sample_every, int enableDriftCorrection,
                                      double *meanOut, double *stderrOut,
                                      int *nSamplesOut)
{
    double *acc = (double *)malloc(N*DIM * sizeof(double));
    int step;
    int halfway = N_steps / 2;
    double sum = 0.0, sumsq = 0.0;
    int nSamples = 0;

    computeAcceleration(vel, acc);
    for (step = 0; step < N_steps; step++) {
        if (enableDriftCorrection && maxBarLengthError() > POS_CORR_TOL) {
            applyPositionCorrection(vel);
            computeAcceleration(vel, acc);
        }
        rattleStep(vel, acc, dt);

        if (step >= halfway && ((step - halfway) % sample_every == 0)) {
            double s = instantaneousStress(vel);
            sum   += s;
            sumsq += s * s;
            nSamples++;
        }
    }
    free(acc);

    double mean = (nSamples > 0) ? sum / nSamples : NAN;
    double var  = 0.0;
    if (nSamples > 1) {
        var = sumsq / nSamples - mean * mean;
        if (var < 0.0) var = 0.0;   /* guard tiny negative from roundoff */
    }
    double se = (nSamples > 1) ? sqrt(var) / sqrt((double)(nSamples - 1)) : 0.0;

    *meanOut     = mean;
    *stderrOut   = se;
    *nSamplesOut = nSamples;
}

/* Driver for the nonlinear (curvature-based) shear modulus measurement,
 * using the single persistent physical velocity `physVel` throughout.
 * physVel is drawn once here, a no-op if an earlier NVE call in this run
 * already initialized it, and never redrawn.
 *
 * The velocity is thermalized first: initializePhysicalVelocity() draws
 * it, then advanceDynamicsToMeasurementState() runs it forward under
 * RATTLE dynamics for N_dyn_steps to reach a genuine dynamical
 * measurement state (R,V). Gamma0 = Gamma(R,V) at the unstrained state is
 * printed for reference; it's a curvature stress, not a modulus, since
 * the modulus itself needs the two-state difference below. For each
 * dgamma in dgammaListStr (comma-separated), computeShearModulusCurvature()
 * gives the genuine two-state finite-difference G_fd -- shear,
 * Lerner-relax to convergence, re-project the same V, difference the
 * stress -- plus the number of Newton relaxation iterations used and the
 * residual bar-length error right after relaxation. One row per dgamma
 * is appended to G_curvature_vs_z.dat, file-locked so it's safe to
 * launch many instances of this binary concurrently for a z-sweep;
 * averaging over realizations into G(z) happens downstream over the
 * accumulated rows.
 *
 * Each dgamma trial starts from the same base state, the post-dynamics
 * measurement state (R,V), not the freshly-drawn velocity from the very
 * first step above. The snapshot taken right after
 * advanceDynamicsToMeasurementState(), and before the sweep, is restored
 * between trials, and once more after the whole sweep. This is pure
 * bookkeeping, a memcpy of saved arrays, and doesn't draw any new random
 * velocity, so it lets several dgamma values be compared apples-to-apples
 * and leaves the network at the measurement state once the sweep
 * finishes.
 */
void runGmod(const char *dgammaListStr, double temperature,
             double dt, int N_dyn_steps)
{
    /* single persistent physical velocity, drawn once, ever */
    initializePhysicalVelocity(temperature);

    /* Run RATTLE dynamics first so the state fed into the curvature
     measurement below is the actual dynamically-evolved constrained
     state (R,V), not the raw just-drawn velocity. */
    if (N_dyn_steps > 0) {
        printf("\n=== Advancing RATTLE dynamics to reach measurement state ===\n");
        printf("  T=%.4g  dt=%.4g  dyn_steps=%d\n", temperature, dt, N_dyn_steps);
        advanceDynamicsToMeasurementState(physVel, dt, N_dyn_steps);
        printf("  Reached (R,V): max_bar_err=%.3e  max_vel_constraint=%.3e\n",
               maxBarLengthError(), maxVelocityConstraint(physVel));
    }

    double z = 2.0 * contacts / DOUBLE_N;
    double Gamma0 = computeCurvatureStress(physVel);   /* Gamma(R,V) at gamma=0, V now dynamically evolved */

    printf("\n=== NONLINEAR (CURVATURE-BASED) SHEAR MODULUS ===\n");
    printf("  z = %.6g   Nc = %d\n", z, contacts);
    printf("  Gamma(R,V) at gamma=0 = %.10g\n", Gamma0);
    printf("  %-12s %-18s %-10s %-14s\n",
           "dgamma", "G_fd", "iters", "residual_err");

    FILE *gout = fopen("G_curvature_vs_z.dat", "a+");
    int haveLock = 0;
    if (gout) {
        int fd = fileno(gout);
        haveLock = (flock(fd, LOCK_EX) == 0);
        if (haveLock) {
            fseek(gout, 0, SEEK_END);
            if (ftell(gout) == 0)
                fprintf(gout, "z\tNc\tdgamma\tG_fd\tGamma0\titers\tresidual_err\n");
        }
    }

    /* Snapshot the base state (geometry + the single physical V) so
     every dgamma trial below starts from the same point, and so the
     network is left exactly as found once the sweep is done */
    double *rx0    = (double *)malloc(N * sizeof(double));
    double *ry0    = (double *)malloc(N * sizeof(double));
    double *vel0   = (double *)malloc(N*DIM * sizeof(double));
    double strain0 = strain, leD0 = leDisplacement;
    memcpy(rx0,  rx,      N*sizeof(double));
    memcpy(ry0,  ry,      N*sizeof(double));
    memcpy(vel0, physVel, N*DIM*sizeof(double));

    char *listCopy = strdup(dgammaListStr);
    char *tok = strtok(listCopy, ",");
    while (tok) {
        double dgamma = atof(tok);
        if (dgamma > 0.0) {
            int iters;
            double residual;
            double G_fd = computeShearModulusCurvature(physVel, dgamma,
                                                         &iters, &residual);
            printf("  %-12.3e %-18.10g %-10d %-14.3e\n",
                   dgamma, G_fd, iters, residual);

            if (haveLock)
                fprintf(gout, "%.10g\t%d\t%.10g\t%.10g\t%.10g\t%d\t%.3e\n",
                        z, contacts, dgamma, G_fd, Gamma0, iters, residual);

            /* restore the base snapshot: next dgamma trial (and the
             caller, once the sweep ends) sees the same state --
             no new random velocity is drawn to do this */
            memcpy(rx,      rx0,  N*sizeof(double));
            memcpy(ry,      ry0,  N*sizeof(double));
            memcpy(physVel, vel0, N*DIM*sizeof(double));
            strain = strain0; leDisplacement = leD0;
            updateGeometry();
        } else {
            fprintf(stderr, "WARNING: skipping non-positive dgamma '%s'\n", tok);
        }
        tok = strtok(NULL, ",");
    }
    free(listCopy);
    free(rx0); free(ry0); free(vel0);

    if (gout) {
        if (haveLock) { fflush(gout); flock(fileno(gout), LOCK_UN); }
        fclose(gout);
    }
}

/* Section 7b: Strain-Ramp Correction Test
 *
 * This section answers a specific empirical question raised by the PI:
 * after imposing a finite strain increment on the bonds, is a *single*
 * geometric correction pass enough to bring the network back onto the
 * constraint manifold (bar error below POS_CORR_TOL), or are several
 * needed? relaxToConstraintManifold() already answers this internally,
 * but it hides the answer behind a hardcoded iteration cap
 * (RELAX_MAX_ITERS) and a very tight tolerance (RELAX_TOL = 1e-12) --
 * there is no record of how many passes were actually used, and no way
 * to change that cap without recompiling.
 *
 * The protocol implemented here, one "strain step" at a time:
 *   1. Impose a finite affine shear increment dgamma (applyAffineShearStep,
 *      already implemented -- unchanged).
 *   2. Apply geometric correction passes ONE AT A TIME
 *      (applyOneGeometricCorrection), checking maxBarLengthError() after
 *      each pass, and stopping as soon as it drops below corr_tol -- or
 *      after max_corrections passes, whichever comes first. Every
 *      pass's residual error is logged, so it's possible to see directly
 *      whether pass #1 alone was already enough.
 *   3. Re-project the persistent physical velocity once onto the
 *      corrected manifold (projectVelocities), the same one-time
 *      pattern already used in computeShearModulusCurvature().
 *   4. Let the system develop dynamically for dyn_steps of RATTLE
 *      (advanceDynamicsToMeasurementState, already implemented --
 *      unchanged; this is where "let the system run long enough to
 *      measure the stress" happens, with its own existing per-step
 *      POS_CORR_TOL drift correction untouched).
 *   5. Measure and log the stress (and bar error) at the end of that
 *      dynamical run.
 * Then repeat for the next strain increment, so the accumulated strain
 * -- and, empirically, the stress -- grows step by step.
 */

/* One Newton/Gauss-Newton correction pass toward the constraint
 * manifold, positions only -- velocity is untouched. This is exactly
 * the loop body used inside relaxToConstraintManifold() and
 * applyPositionCorrection(), pulled out as its own function so a caller
 * can apply one pass at a time and decide whether another is needed,
 * instead of having the iteration count and stopping tolerance baked
 * in. Solves N f = dr for the current bar-length errors dr, then shifts
 * positions by dX = -S^T f (minimum-norm correction) and rebuilds nij[]. */
void applyOneGeometricCorrection(void)
{
    double *rCur = (double *)malloc(contacts * sizeof(double));
    double *dr   = (double *)malloc(contacts * sizeof(double));
    double *dX   = (double *)malloc(N*DIM    * sizeof(double));
    int l, i;

    computeBarLengths(rCur);
    for (l = 0; l < contacts; l++)
        dr[l] = rCur[l] - restLength[l];
    free(rCur);

    double *f = solveN(dr);
    esTransposeOnVector(f, dX);
    free(f); free(dr);

    for (i = 0; i < N; i++) {
        double newrx = rx[i] - dX[DIM*i + X_COMP] * invL;
        double newry = ry[i] - dX[DIM*i + Y_COMP] * invL;
        if      (newrx >= LENGTH) newrx -= LENGTH;
        else if (newrx <  0.0)   newrx += LENGTH;
        if      (newry >= LENGTH) newry -= LENGTH;
        else if (newry <  0.0)   newry += LENGTH;
        rx[i] = newrx; ry[i] = newry;
    }
    free(dX);

    updateGeometry();
}

/* Strain-ramp driver. Runs n_strain_steps repetitions of the protocol
 * described above, appending one row per correction pass to
 * corrTraceFileName (so the convergence of the correction loop itself
 * can be inspected) and one row per strain step to summaryFileName (so
 * stress vs. accumulated strain can be inspected). Both files are
 * opened fresh ("w") at the start of the run -- this is a diagnostic
 * test driver, not meant to be run concurrently across many processes
 * the way runGmod() is, so no file locking / append-and-share logic is
 * needed here. */
void runStrainTest(int n_strain_steps, double dgamma, int enableCorrection,
                    int max_corrections, double corr_tol, int dyn_steps,
                    int sample_every, double temperature, double dt,
                    char *summaryFileName, char *corrTraceFileName)
{
    initializePhysicalVelocity(temperature);

    FILE *fsum  = fopen(summaryFileName, "w");
    FILE *ftrace = fopen(corrTraceFileName, "w");
    if (!fsum || !ftrace) {
        printf("couldn't open output files for straintest\n");
        if (fsum) fclose(fsum);
        if (ftrace) fclose(ftrace);
        return;
    }

    fprintf(fsum, "# rigid-bar strain-ramp correction test  N=%d  z0=%.4g\n",
            N, 2.0*contacts/DOUBLE_N);
    fprintf(fsum, "# enableCorrection=%d (0 = no correction ANYWHERE -- neither\n"
                  "#   the pre-dynamics correction loop nor the internal\n"
                  "#   per-step drift check inside the dynamics window are run)\n",
            enableCorrection);
    fprintf(fsum, "# dgamma=%.6g  max_corrections=%d  corr_tol=%.3e  "
                  "dyn_steps=%d  sample_every=%d  T=%.4g  dt=%.4g\n",
            dgamma, max_corrections, corr_tol, dyn_steps, sample_every,
            temperature, dt);
    fprintf(fsum, "# col: strain_step  cumulative_strain  n_corrections_used  "
                  "converged  bar_err_after_correction  sigma_xy_preCorrection  "
                  "sigma_xy_postCorrection  sigma_xy_endpoint  sigma_xy_mean  "
                  "sigma_xy_stderr  n_samples  pressure_postdyn  bar_err_postdyn  "
                  "vel_constraint_postdyn  KE_postdyn\n");

    fprintf(ftrace, "# strain_step  cumulative_strain  corr_iter  max_bar_err\n");

    printf("\n=== Strain-ramp correction test ===\n");
    printf("  dgamma=%.4g  max_corrections=%d  corr_tol=%.3e  dyn_steps=%d  "
           "sample_every=%d\n",
           dgamma, max_corrections, corr_tol, dyn_steps, sample_every);

    int k;
    for (k = 1; k <= n_strain_steps; k++) {

        /* 1. Impose the strain increment */
        applyAffineShearStep(dgamma);

        /* Gamma measured on the just-sheared geometry, BEFORE any
         * correction is applied -- same physVel that survived from the
         * end of the previous strain step's dynamics run, on the
         * geometry immediately after the affine map (updateGeometry()
         * inside applyAffineShearStep() already refreshed nij[], so this
         * is well-defined even though the bars are momentarily off their
         * rest length). This isolates whether the correction step itself
         * changes the measured stress, holding velocity fixed. */
        double sigma_preCorrection = instantaneousStress(physVel);

        /* 2. Correction passes, one at a time, logged individually.
         * Skipped entirely when enableCorrection=0 -- this is the
         * "watch it happen with no correction anywhere" mode. */
        double barErr = maxBarLengthError();
        fprintf(ftrace, "%d\t%.10g\t%d\t%.6e\n", k, strain, 0, barErr);

        int n_corr = 0;
        if (enableCorrection) {
            while (barErr > corr_tol && n_corr < max_corrections) {
                applyOneGeometricCorrection();
                n_corr++;
                barErr = maxBarLengthError();
                fprintf(ftrace, "%d\t%.10g\t%d\t%.6e\n", k, strain, n_corr, barErr);
            }
        }
        int converged = (barErr <= corr_tol) ? 1 : 0;

        /* 3. Re-project the persistent physical velocity once onto the
         * corrected manifold -- same one-shot pattern already used in
         * computeShearModulusCurvature(). Also skipped when
         * enableCorrection=0, so velocity is left exactly as RATTLE's
         * own dynamics produced it, with no re-projection step at all.
         * Gamma measured again right after, same velocity field (now
         * reprojected iff enableCorrection), corrected-or-not geometry
         * -- the direct counterpart to sigma_preCorrection above. */
        if (enableCorrection)
            projectVelocities(physVel);
        double sigma_postCorrection = instantaneousStress(physVel);

        /* 4. Let the system develop dynamically, sampling Gamma over
         * the second half of the window and averaging, so a single
         * noisy endpoint snapshot can be compared directly against a
         * time-averaged estimate. The internal per-step drift check is
         * likewise gated by enableCorrection -- with it off, this is a
         * genuinely uncorrected RATTLE run, not just a run with the
         * pre-dynamics loop skipped (that alone still let the internal
         * drift check silently patch things up mid-window). */
        double sigma_mean = NAN, sigma_stderr = 0.0;
        int n_samples = 0;
        if (dyn_steps > 0) {
            advanceDynamicsAndAverageStress(physVel, dt, dyn_steps, sample_every,
                                             enableCorrection,
                                             &sigma_mean, &sigma_stderr, &n_samples);
        } else {
            sigma_mean = instantaneousStress(physVel);
        }

        /* 5. Measure */
        double sigma_endpoint = instantaneousStress(physVel);
        double pressure_postdyn = computePressure(physVel);
        double barErr_postdyn = maxBarLengthError();
        double velCon_postdyn = maxVelocityConstraint(physVel);
        double KE_postdyn = kineticEnergy(physVel);

        fprintf(fsum, "%d\t%.10g\t%d\t%d\t%.6e\t%.10g\t%.10g\t%.10g\t%.10g\t%.6e\t%d\t"
                      "%.10g\t%.6e\t%.6e\t%.10g\n",
                k, strain, n_corr, converged, barErr, sigma_preCorrection,
                sigma_postCorrection, sigma_endpoint, sigma_mean, sigma_stderr,
                n_samples, pressure_postdyn, barErr_postdyn, velCon_postdyn,
                KE_postdyn);
        fflush(fsum); fflush(ftrace);

        printf("  strain_step %3d  cum_strain=%.4g  n_corr=%2d  converged=%d  "
               "barErr=%.2e  sigma(preCorr)=%.4g  sigma(postCorr)=%.4g  "
               "sigma(mean)=%.4g +/- %.2e (n=%d)  barErr(postdyn)=%.2e\n",
               k, strain, n_corr, converged, barErr, sigma_preCorrection,
               sigma_postCorrection, sigma_mean, sigma_stderr,
               n_samples, barErr_postdyn);
    }

    fclose(fsum);
    fclose(ftrace);
    printf("\nStrain-ramp test complete.\n  summary  -> %s\n  corr trace -> %s\n",
           summaryFileName, corrTraceFileName);
}

/* Section 8: Main */
int main(int argc, char **argv)
{
    char networkFile[1024], outFile[1024], *p;
    double T = 0.01, dt = 0.002;
    int N_steps = 200000, out_every = 1000;
    int doNVE = 0, doGmod = 0, doStrainTest = 0;
    char gmodList[1024] = "1e-3,1e-4,1e-5,1e-6";
    double Tgmod = 0.01;
    double dtGmod = 0.002;
    int dynStepsGmod = 10000;

    /* straintest defaults -- deliberately small so a quick interactive
     * run stays well under a minute even on a laptop; override from the
     * command line for a real study. */
    int    nStrainSteps    = 5;
    double dgammaStrain    = 1e-3;
    int    enableCorrection = 1;   /* 0 = no correction anywhere (neither the
                                       pre-dynamics loop nor the per-step
                                       drift check during dynamics) -- lets
                                       you watch the raw, uncorrected drift */
    int    maxCorrections  = 50;
    double corrTol         = POS_CORR_TOL;   /* 1e-8, same tolerance runNVE's
                                                 drift check already uses */
    int    dynStepsStrain  = 50;
    int    sampleEvery     = 5;   /* stride, in RATTLE steps, between Gamma
                                     samples taken over the second half of
                                     the dyn_steps window */
    double Tstrain         = 0.01;
    double dtStrain        = 0.002;

    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s <network.dat> nve [T] [dt] [steps] [out_every]\n"
            "       %s <network.dat> gmod [dgamma_list] [T] [dt] [dyn_steps]\n"
            "       %s <network.dat> nvegmod [dgamma_list] [T] [dt] [steps] [out_every]\n"
            "                 (logs the RATTLE run to <network>_nve.dat, same as nve,\n"
            "                  then measures G(z) from the exact state that run reached\n"
            "                  -- no separate, unlogged thermalization pass)\n"
            "       %s <network.dat> straintest [n_strain_steps] [dgamma] [enable_correction] [max_corrections] [corr_tol] [dyn_steps] [sample_every] [T] [dt]\n"
            "                 (enable_correction=0: no correction ANYWHERE -- neither the\n"
            "                  pre-dynamics correction loop nor the internal per-step drift\n"
            "                  check during the dynamics window run at all, so you see the\n"
            "                  raw, uncorrected drift in bar error and stress. Use this first\n"
            "                  to see the phenomenon with no correction involved, then rerun\n"
            "                  with enable_correction=1 (default) to see how many corrections\n"
            "                  it takes to stay under corr_tol.)\n"
            "                 (repeatedly: impose dgamma, apply up to max_corrections\n"
            "                  single-pass geometric corrections logging each pass's\n"
            "                  residual bar error, let the system run dyn_steps of\n"
            "                  RATTLE while sampling Gamma every sample_every steps over\n"
            "                  the second half of that window (discarding the first half\n"
            "                  as a transient) to get a mean +/- standard error, in\n"
            "                  addition to the raw endpoint snapshot -- writes\n"
            "                  <network>_straintest_summary.dat and\n"
            "                  <network>_straintest_corrtrace.dat)\n",
            argv[0], argv[0], argv[0], argv[0]);
        return 1;
    }

    /* Seed the RNG once per process. time(NULL) alone repeats if several instances of this binary launch within the same second
    (e.g. a Python loop over many network files), so mix in the PID to keep independent runs independent. */
    srand((unsigned int)(time(NULL) ^ (unsigned int)getpid()));

    strncpy(networkFile, argv[1], sizeof(networkFile)-1);

    if (argc >= 3 && strcmp(argv[2], "nve") == 0) {
        doNVE = 1;
        if (argc >= 4) T        = atof(argv[3]);
        if (argc >= 5) dt       = atof(argv[4]);
        if (argc >= 6) N_steps  = atoi(argv[5]);
        if (argc >= 7) out_every = atoi(argv[6]);
    } else if (argc >= 3 && strcmp(argv[2], "gmod") == 0) {
        doGmod = 1;
        if (argc >= 4) strncpy(gmodList, argv[3], sizeof(gmodList)-1);
        if (argc >= 5) Tgmod        = atof(argv[4]);
        if (argc >= 6) dtGmod       = atof(argv[5]);
        if (argc >= 7) dynStepsGmod = atoi(argv[6]);
    } else if (argc >= 3 && strcmp(argv[2], "nvegmod") == 0) {
        doNVE = 1;
        doGmod = 1;
        if (argc >= 4) strncpy(gmodList, argv[3], sizeof(gmodList)-1);
        if (argc >= 5) T         = atof(argv[4]);
        if (argc >= 6) dt        = atof(argv[5]);
        if (argc >= 7) N_steps   = atoi(argv[6]);
        if (argc >= 8) out_every = atoi(argv[7]);
        /* Same T, dt as the logged NVE run above -- this is one dynamical
         * trajectory, not two. dynStepsGmod=0 makes runGmod() skip its own
         * advanceDynamicsToMeasurementState() call entirely (see the
         * N_dyn_steps>0 guard in runGmod()), so it measures G(z) from
         * exactly the state runNVE() just finished logging, instead of
         * quietly running more unlogged dynamics on top of it. */
        Tgmod = T;
        dtGmod = dt;
        dynStepsGmod = 0;
    } else if (argc >= 3 && strcmp(argv[2], "straintest") == 0) {
        doStrainTest = 1;
        if (argc >=  4) nStrainSteps    = atoi(argv[3]);
        if (argc >=  5) dgammaStrain    = atof(argv[4]);
        if (argc >=  6) enableCorrection = atoi(argv[5]);
        if (argc >=  7) maxCorrections  = atoi(argv[6]);
        if (argc >=  8) corrTol         = atof(argv[7]);
        if (argc >=  9) dynStepsStrain  = atoi(argv[8]);
        if (argc >= 10) sampleEvery     = atoi(argv[9]);
        if (argc >= 11) Tstrain         = atof(argv[10]);
        if (argc >= 12) dtStrain        = atof(argv[11]);
        if (sampleEvery < 1) sampleEvery = 1;
    }

    /* read network */
    if (!readNetwork(networkFile)) return 1;

    /* build output filename by stripping .dat from input and appending _nve.dat */
    snprintf(outFile, sizeof(outFile), "%s", networkFile);
    p = strstr(outFile, ".dat");
    if (p) *p = '\0';
    strncat(outFile, "_nve.dat", sizeof(outFile) - strlen(outFile) - 1);

    /* NVE dynamics */
    if (doNVE) {
        verboseCG = 0;
        runNVE(T, dt, N_steps, out_every, outFile);
    }

    /* nonlinear (curvature-based) G(z) */
    if (doGmod) {
        verboseCG = 0;
        if (doNVE)
            printf("\n(gmod measuring the state just logged above -- "
                   "dyn_steps=0, no additional dynamics)\n");
        runGmod(gmodList, Tgmod, dtGmod, dynStepsGmod);
    }

    /* strain-ramp correction test */
    if (doStrainTest) {
        char summaryFile[1024], traceFile[1024];
        verboseCG = 0;
        snprintf(summaryFile, sizeof(summaryFile), "%s", networkFile);
        p = strstr(summaryFile, ".dat");
        if (p) *p = '\0';
        snprintf(traceFile, sizeof(traceFile), "%s", summaryFile);
        strncat(summaryFile, "_straintest_summary.dat",
                sizeof(summaryFile) - strlen(summaryFile) - 1);
        strncat(traceFile, "_straintest_corrtrace.dat",
                sizeof(traceFile) - strlen(traceFile) - 1);
        runStrainTest(nStrainSteps, dgammaStrain, enableCorrection, maxCorrections,
                      corrTol, dynStepsStrain, sampleEvery, Tstrain, dtStrain,
                      summaryFile, traceFile);
    }

    /* physVel (if it was ever initialized, by runNVE and/or runGmod) is
     * the single persistent physical velocity for the whole run; free
     * it exactly once, here, at the very end. */
    if (physVel) { free(physVel); physVel = NULL; }

    freeNetwork();
    return 0;
}

/* _fftw_dst_shim.c
 *
 * Thin C shim exposing exactly the FFTW3 calls src/EllipticSolver2d.cc uses
 * for its 2D discrete sine transform (DST-I / FFTW_RODFT00), so the Python
 * port (py/elliptic_solver_2d.py) can call the SAME FFTW3 library, with the
 * SAME planning flag (FFTW_EXHAUSTIVE), as the C++ reference code -- not an
 * equivalent-result substitute (e.g. scipy.fft), but literally the same
 * planner searching the same list of candidate sine-transform algorithms
 * ("codelets") for the fastest one on this machine, for this problem size.
 *
 * Mirrors, line for line, the C++ constructor/sinTransform in
 * src/EllipticSolver2d.cc:
 *     _FFTWPlan = fftw_plan_r2r_2d( nx-1, ny-1, _fft, _fft,
 *         FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE );
 *     ...
 *     fftw_execute( _FFTWPlan );
 *
 * Built on demand by py/_fftw_native.py (via clang, into build/, which is
 * gitignored -- this .c file is the only thing committed).
 */

#include <fftw3.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n0, n1;
    fftw_plan plan;
    double *buf;
} dst2d_plan_t;

/* Build a plan for an n0 x n1 real-to-real DST-I (FFTW_RODFT00) transform,
 * using FFTW_EXHAUSTIVE -- i.e. actually timing every algorithm FFTW knows
 * for this transform type/size and keeping whichever measured fastest.
 * Matches EllipticSolver2d's constructor exactly (same call, same flag,
 * same transform kind, same array size convention: n0=nx-1, n1=ny-1). */
void *dst2d_create(int n0, int n1) {
    double *buf = (double *)fftw_malloc(sizeof(double) * (size_t)n0 * (size_t)n1);
    fftw_plan p = fftw_plan_r2r_2d(n0, n1, buf, buf, FFTW_RODFT00, FFTW_RODFT00,
                                    FFTW_EXHAUSTIVE);
    dst2d_plan_t *h = (dst2d_plan_t *)malloc(sizeof(dst2d_plan_t));
    h->n0 = n0;
    h->n1 = n1;
    h->plan = p;
    h->buf = buf;
    return (void *)h;
}

/* Run the plan built by dst2d_create on `data` (n0*n1 doubles, row-major,
 * in place) -- matches EllipticSolver2d::sinTransform's copy-in /
 * fftw_execute / copy-out sequence. */
void dst2d_execute(void *handle, double *data) {
    dst2d_plan_t *h = (dst2d_plan_t *)handle;
    size_t nbytes = sizeof(double) * (size_t)h->n0 * (size_t)h->n1;
    memcpy(h->buf, data, nbytes);
    fftw_execute(h->plan);
    memcpy(data, h->buf, nbytes);
}

/* Release a plan created by dst2d_create -- matches
 * EllipticSolver2d::~EllipticSolver2d's fftw_destroy_plan call. */
void dst2d_destroy(void *handle) {
    dst2d_plan_t *h = (dst2d_plan_t *)handle;
    fftw_destroy_plan(h->plan);
    fftw_free(h->buf);
    free(h);
}

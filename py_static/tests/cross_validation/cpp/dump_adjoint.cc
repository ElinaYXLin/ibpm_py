// dump_adjoint.cc
//
// Raw-value dump driver for cross-validating py/vector_operations.py's
// FluxToYVelocity / YVelocityToFlux against the C++ reference. Reproduces
// the adjointness check from test/CheckAdjoint.cc (verify that
// YVelocityToFlux is the adjoint of FluxToYVelocity, computed brute-force
// per basis vector, across a 3-level multigrid), but prints "mag" and
// "err2" as one full-precision value per line instead of via the original's
// printY() (default 6-significant-digit cout formatting).
//
// Build (from repo root):
//   g++ -Wall -O2 -I src -c py/tests/cross_validation/cpp/dump_adjoint.cc -o /tmp/dump_adjoint.o
//   g++ /tmp/dump_adjoint.o -L build -libpm -lfftw3 -lm -o /tmp/dump_adjoint
// Run:
//   /tmp/dump_adjoint > /tmp/dump_adjoint.out

#include <iostream>
#include <iomanip>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

// Compute the adjoint of YVelocityToFlux, brute force:
// x(lev,i,j) = < q, p(lev,i,j) > / <e,e> where e is an orthogonal basis function
void AdjYVelToFlux( const Flux& q, Scalar& x ) {
    Scalar e( x.getGrid() );
    x = 0.;
    for (int lev=0; lev<x.Ngrid(); ++lev) {
        for (int i=1; i<x.Nx(); ++i) {
            for (int j=1; j<x.Ny(); ++j) {
                e = 0;
                e(lev,i,j) = 1.;
                Flux p(q.getGrid());
                p = 0;
                YVelocityToFlux( e, p );
                double a = InnerProduct( q, p );
                double normsq = InnerProduct( e, e );
                x(lev,i,j) = a / normsq;
            }
        }
    }
}

void dumpFluxY(const char* name, const Flux& q) {
    cout << name << endl;
    for (int lev=0; lev<q.Ngrid(); ++lev)
        for (int i=0; i<q.Nx(); ++i)
            for (int j=0; j<=q.Ny(); ++j)
                cout << setprecision(17) << q(lev,Y,i,j) << endl;
}

int main() {
    int nx=8, ny=8, ngrid=3;
    double length=0.1;
    Grid grid( nx, ny, ngrid, length, 0., 0. );
    Flux mag(grid);
    Flux err2(grid);
    mag = 0.;
    err2 = 0.;

    // Loop over all fluxes in the Y direction
    for (int lev=0; lev<ngrid; ++lev) {
        for (int ind=err2.begin(Y); ind < err2.end(Y); ++ind) {
            Flux e(grid);
            e = 0;
            e(lev,ind) = 1;
            Scalar x(grid);
            Scalar x0(grid);
            x = 0;
            x0 = 0;
            FluxToYVelocity( e, x0 );
            AdjYVelToFlux( e, x );
            mag(lev,ind) = InnerProduct(x0,x0);
            err2(lev,ind) = InnerProduct(x-x0,x-x0);
        }
    }
    dumpFluxY("DUMP_mag", mag);
    dumpFluxY("DUMP_err2", err2);
    return 0;
}

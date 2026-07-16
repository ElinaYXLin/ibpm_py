// dump_fluxtox.cc
//
// Raw-value dump driver for cross-validating py/vector_operations.py's
// FluxToXVelocity / XVelocityToFlux against the C++ reference. Derived from
// test/CheckFluxToX.cc, with the interactive `cin.get(c)` debug loop removed
// (it isn't a real automated check -- it printed each grid point's
// dependency coefficients and waited for a keypress) and replaced by a
// fixed spot-check of computeDependencies() at 7 representative (lev,i,j)
// points spanning the 3 multigrid levels and interior/border/corner regions.
//
// Build (from repo root):
//   g++ -Wall -O2 -I src -c py/tests/cross_validation/cpp/dump_fluxtox.cc -o /tmp/dump_fluxtox.o
//   g++ /tmp/dump_fluxtox.o -L build -libpm -lfftw3 -lm -o /tmp/dump_fluxtox
// Run:
//   /tmp/dump_fluxtox > /tmp/dump_fluxtox.out

#include <iostream>
#include <iomanip>
#include "ibpm.h"
using namespace std;
using namespace ibpm;

void dumpScalar(const char* name, const Scalar& s) {
    cout << name << endl;
    for (int lev=0; lev<s.Ngrid(); ++lev)
        for (int i=1; i<s.Nx(); ++i)
            for (int j=1; j<s.Ny(); ++j)
                cout << setprecision(17) << s(lev,i,j) << endl;
}
void dumpFluxX(const char* name, const Flux& q) {
    cout << name << endl;
    for (int lev=0; lev<q.Ngrid(); ++lev)
        for (int i=0; i<=q.Nx(); ++i)
            for (int j=0; j<q.Ny(); ++j)
                cout << setprecision(17) << q(lev,X,i,j) << endl;
}

// In FluxToXVelocity, scalar x(lev,i,j) depends on nearby flux values.
// Set q to the coefficient of each flux value that influences x(lev,i,j).
void computeDependencies(int lev, int i, int j, Flux& q) {
    q = 0.;
    for (int l=0; l < q.Ngrid(); ++l) {
        for (int ind=q.begin(X); ind != q.end(X); ++ind) {
            Flux e( q.getGrid() );
            Scalar x( q.getGrid() );
            e = 0.;
            e(l,ind) = 1.;
            FluxToXVelocity(e, x);
            q(l,ind) = x(lev,i,j);
        }
    }
}

int main() {
    int nx=8, ny=8, ngrid=3;
    double length=8;
    Grid grid( nx, ny, ngrid, length, 0., 0. );

    // Round-trip: constant field
    Scalar u(grid);
    u = 1.;
    Flux p(grid);
    XVelocityToFlux( u, p );
    FluxToXVelocity( p, u );
    dumpFluxX("DUMP_ToFlux_u", p);
    dumpScalar("DUMP_roundtrip_u", u);

    // Spot-check computeDependencies at a handful of (lev,i,j)
    int pts[][3] = { {0,1,1}, {0,4,4}, {0,7,7}, {1,2,3}, {1,4,4}, {2,1,1}, {2,4,4} };
    for (auto& pt : pts) {
        Flux q(grid);
        computeDependencies(pt[0], pt[1], pt[2], q);
        cout << "DUMP_dep_" << pt[0] << "_" << pt[1] << "_" << pt[2] << endl;
        for (int lev=0; lev<q.Ngrid(); ++lev)
            for (int i=0; i<=q.Nx(); ++i)
                for (int j=0; j<q.Ny(); ++j)
                    cout << setprecision(17) << q(lev,X,i,j) << endl;
    }
    return 0;
}

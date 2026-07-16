// dump_curl.cc
//
// Raw-value dump driver for cross-validating py/vector_operations.py's Curl
// against the C++ reference implementation. Reproduces the exact scenario
// from test/CheckCurl.cc (same grid, same analytic fields), but prints one
// full-precision value per line under a "DUMP_<name>" marker instead of
// Scalar::print()/Flux::print()'s fixed-width grid layout, so the output is
// trivially diffable against python/dump_curl.py's output via compare_dumps.py.
//
// Build (from repo root):
//   g++ -Wall -O2 -I src -c py/tests/cross_validation/cpp/dump_curl.cc -o /tmp/dump_curl.o
//   g++ /tmp/dump_curl.o -L build -libpm -lfftw3 -lm -o /tmp/dump_curl
// Run:
//   /tmp/dump_curl > /tmp/dump_curl.out

#include <iostream>
#include <iomanip>
#include "ibpm.h"

using namespace std;
using namespace ibpm;

void AdjointOfScalarCurl( const Flux& q, Scalar& f ) {
    for (int lev=0; lev<q.Ngrid(); ++lev) {
        for (int i=1; i<q.Nx(); ++i) {
            for (int j=1; j<q.Ny(); ++j) {
                Scalar e(q.getGrid());
                e = 0;
                e(lev,i,j) = 1.;
                Flux curlE(q.getGrid()) ;
                Curl( e, curlE );
                f(lev,i,j) = InnerProduct( q, curlE ) / InnerProduct( e, e );
            }
        }
    }
}

void dumpScalar(const char* name, const Scalar& s) {
    cout << name << endl;
    for (int lev=0; lev<s.Ngrid(); ++lev)
        for (int i=1; i<s.Nx(); ++i)
            for (int j=1; j<s.Ny(); ++j)
                cout << setprecision(17) << s(lev,i,j) << endl;
}
void dumpFlux(const char* name, const Flux& q) {
    cout << name << endl;
    for (int lev=0; lev<q.Ngrid(); ++lev) {
        for (int i=0; i<=q.Nx(); ++i)
            for (int j=0; j<q.Ny(); ++j)
                cout << setprecision(17) << q(lev,X,i,j) << endl;
        for (int i=0; i<q.Nx(); ++i)
            for (int j=0; j<=q.Ny(); ++j)
                cout << setprecision(17) << q(lev,Y,i,j) << endl;
    }
}

int main() {
    int nx=8, ny=8, ngrid=2;
    double length=8;
    Grid grid( nx, ny, ngrid, length, -length/2, -length*ny/nx/2 );
    Scalar u(grid);
    Scalar v(grid);
    Flux q(grid);

    for (int lev=0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                u(lev,i,j) = grid.getXEdge(lev,i);
                v(lev,i,j) = grid.getYEdge(lev,j);
            }
        }
    }

    for (int lev=0; lev<ngrid; ++lev) {
        double dx = grid.Dx(lev);
        for (int i=0; i<=nx; ++i)
            for (int j=0; j<ny; ++j)
                q(lev,X,i,j) = -2*q.y(lev,X,j)*q.y(lev,X,j) * dx;
        for (int i=0; i<nx; ++i)
            for (int j=0; j<=ny; ++j)
                q(lev,Y,i,j) = q.x(lev,Y,i)*q.x(lev,Y,i) * dx;
    }

    Flux curlU(grid);
    Flux curlV(grid);
    Curl( u, curlU );
    Curl( v, curlV );

    Scalar omega(grid);
    Curl( q, omega );
    Scalar omega2(grid);
    AdjointOfScalarCurl( q, omega2 );
    Scalar err = omega - omega2;

    dumpScalar("DUMP_u", u);
    dumpFlux("DUMP_curlU", curlU);
    dumpScalar("DUMP_v", v);
    dumpFlux("DUMP_curlV", curlV);
    dumpFlux("DUMP_q", q);
    dumpScalar("DUMP_omega", omega);
    dumpScalar("DUMP_omega2", omega2);
    dumpScalar("DUMP_err", err);
    return 0;
}

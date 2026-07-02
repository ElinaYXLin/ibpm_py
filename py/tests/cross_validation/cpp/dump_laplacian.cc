// dump_laplacian.cc
//
// Raw-value dump driver for cross-validating py/vector_operations.py's
// Laplacian (and Scalar.coarsify()) against the C++ reference. Reproduces
// the checkLCEqualsCLC scenario from test/CheckLaplacian.cc: for each grid
// point, checks whether L(coarsify(u)) == coarsify(L(coarsify(u))), where u
// is a unit basis vector. This mostly does NOT hold (it's testing a genuine
// mathematical property of the discretization, not asserting correctness)
// -- the point of this cross-validation is that the exact pass/fail PATTERN
// must match between C++ and Python.
//
// Build (from repo root):
//   g++ -Wall -O2 -I src -c py/tests/cross_validation/cpp/dump_laplacian.cc -o /tmp/dump_laplacian.o
//   g++ /tmp/dump_laplacian.o -L build -libpm -lfftw3 -lm -o /tmp/dump_laplacian
// Run:
//   /tmp/dump_laplacian > /tmp/dump_laplacian.out

#include <iostream>
#include <iomanip>
#include "ibpm.h"
using namespace std;
using namespace ibpm;

PoissonSolver* solver;

bool equal( const Scalar& u, const Scalar& v ) {
    static double tol=1e-12;
    double err = 0;
    for (int lev=0; lev<u.Ngrid(); ++lev)
        for (int i=1; i<u.Nx(); ++i)
            for (int j=1; j<u.Ny(); ++j)
                err += fabs( u(lev,i,j) - v(lev,i,j) );
    return err < tol;
}

bool checkLCEqualsCLC( const Grid& grid, int lev, int i, int j ) {
    Scalar u(grid);
    u = 0;
    u(lev,i,j) = 1;
    u.coarsify();
    Scalar LCu(grid);
    Laplacian( u, LCu );
    Scalar CLCu = LCu;
    CLCu.coarsify();
    return equal( CLCu, LCu );
}

int main() {
    int nx=8, ny=8, ngrid=2;
    double length=8;
    Grid grid( nx, ny, ngrid, length, -length/2, -length*ny/nx/2 );
    solver = new PoissonSolver( grid );

    Scalar status(grid);
    status = 2;

    for (int lev = 0; lev<ngrid; ++lev) {
        for (int i=1; i<nx; ++i) {
            for (int j=1; j<ny; ++j) {
                if ( checkLCEqualsCLC( grid, lev, i, j ) ) {
                    status(lev,i,j) = 0;
                } else {
                    status(lev,i,j) = 1;
                }
            }
        }
    }
    cout << "DUMP_status" << endl;
    for (int lev=0; lev<ngrid; ++lev)
        for (int i=1; i<nx; ++i)
            for (int j=1; j<ny; ++j)
                cout << setprecision(17) << status(lev,i,j) << endl;
    return 0;
}

"""
Read and validate IBPM binary restart files.

Binary format (from State.cc save()):
  - Grid info (6 ints + doubles):
    - nx, ny, ngrid (int x 3)
    - dx, x0, y0 (double x 3)
  - Geometry info:
    - numPoints (int)
  - Flux q data: all levels, then all indices
  - Scalar omega data: all levels, interior points only (1..nx-1, 1..ny-1)
  - BoundaryVector f: all boundary points
  - Timestep and time info

Flux storage (from Flux.h):
  - X-fluxes: (nx+1) x ny, indexed as i*(ny) + j
  - Y-fluxes: nx x (ny+1), indexed as i*(ny+1) + j
  - Total fluxes per level: (nx+1)*ny + nx*(ny+1)

Scalar storage (from Scalar.h):
  - Interior points only: (nx-1) x (ny-1), 1-indexed loop i=1..nx-1, j=1..ny-1

BoundaryVector storage (from BoundaryVector.h):
  - numPoints vector pairs: [fx1, fy1, fx2, fy2, ...]
"""

import struct
import numpy as np
from pathlib import Path


class IBPMRestart:
    """Reader for IBPM binary restart files."""
    
    def __init__(self, filename):
        self.filename = filename
        self.data = {}
        self.read()
    
    def read(self):
        """Read binary restart file."""
        with open(self.filename, 'rb') as f:
            # Read Grid info
            nx, ny, ngrid = struct.unpack('3i', f.read(3 * 4))
            dx, x0, y0 = struct.unpack('3d', f.read(3 * 8))
            
            # Read Geometry info
            numPoints = struct.unpack('i', f.read(4))[0]
            
            self.data['nx'] = nx
            self.data['ny'] = ny
            self.data['ngrid'] = ngrid
            self.data['dx'] = dx
            self.data['x0'] = x0
            self.data['y0'] = y0
            self.data['numPoints'] = numPoints
            
            # Calculate flux dimensions
            # X-fluxes: (nx+1) x ny
            # Y-fluxes: nx x (ny+1)
            num_x_fluxes = (nx + 1) * ny
            num_fluxes_per_level = num_x_fluxes + nx * (ny + 1)
            
            # Read Flux q
            q_flat = np.frombuffer(
                f.read(ngrid * num_fluxes_per_level * 8),
                dtype=np.float64
            )
            self.data['q_flat'] = q_flat
            
            # Reshape: (ngrid, num_fluxes_per_level)
            self.data['q'] = q_flat.reshape((ngrid, num_fluxes_per_level))
            
            # Read Scalar omega (interior points only)
            num_interior = (nx - 1) * (ny - 1)
            omega_flat = np.frombuffer(
                f.read(ngrid * num_interior * 8),
                dtype=np.float64
            )
            self.data['omega_flat'] = omega_flat
            
            # Reshape omega with padding for boundaries (all zeros)
            # Storage is just interior, but we'll create full arrays with zeros at boundaries
            self.data['omega'] = np.zeros((ngrid, nx, ny))
            idx = 0
            for lev in range(ngrid):
                for i in range(1, nx - 1):
                    for j in range(1, ny - 1):
                        self.data['omega'][lev, i, j] = omega_flat[idx]
                        idx += 1
            
            # Read BoundaryVector f (forces)
            f_data = np.frombuffer(
                f.read(numPoints * 2 * 8),
                dtype=np.float64
            )
            # Reshape: (numPoints, 2) with [fx, fy]
            self.data['f'] = f_data.reshape((numPoints, 2))
            
            # Read timestep and time
            timestep = struct.unpack('i', f.read(4))[0]
            time = struct.unpack('d', f.read(8))[0]
            
            self.data['timestep'] = timestep
            self.data['time'] = time
    
    def get_flux_xy(self):
        """
        Extract X and Y flux components as separate arrays.
        
        Returns:
            q_x: (ngrid, nx+1, ny) array of x-fluxes
            q_y: (ngrid, nx, ny+1) array of y-fluxes
        """
        nx = self.data['nx']
        ny = self.data['ny']
        ngrid = self.data['ngrid']
        q = self.data['q']
        
        num_x_fluxes = (nx + 1) * ny
        
        q_x = np.zeros((ngrid, nx + 1, ny))
        q_y = np.zeros((ngrid, nx, ny + 1))
        
        for lev in range(ngrid):
            # Extract X-fluxes: indices 0..num_x_fluxes-1
            for i in range(nx + 1):
                for j in range(ny):
                    idx = i * ny + j
                    q_x[lev, i, j] = q[lev, idx]
            
            # Extract Y-fluxes: indices num_x_fluxes..end
            for i in range(nx):
                for j in range(ny + 1):
                    idx = num_x_fluxes + i * (ny + 1) + j
                    q_y[lev, i, j] = q[lev, idx]
        
        return q_x, q_y
    
    def get_omega(self):
        """Get vorticity field (full array with zero boundary)."""
        return self.data['omega']
    
    def get_forces(self):
        """Get boundary forces."""
        return self.data['f']
    
    def info(self):
        """Print summary of restart file."""
        print(f"Restart file: {self.filename}")
        print(f"  Grid: nx={self.data['nx']}, ny={self.data['ny']}, ngrid={self.data['ngrid']}")
        print(f"  Grid spacing: dx={self.data['dx']:.6f}")
        print(f"  Grid origin: x0={self.data['x0']:.6f}, y0={self.data['y0']:.6f}")
        print(f"  Boundary points: {self.data['numPoints']}")
        print(f"  Timestep: {self.data['timestep']}, Time: {self.data['time']:.6f}")
        print(f"  Flux shape: {self.data['q'].shape}")
        print(f"  Omega shape: {self.data['omega'].shape}")
        print(f"  Forces shape: {self.data['f'].shape}")
        
        q_x, q_y = self.get_flux_xy()
        omega = self.get_omega()
        f = self.get_forces()
        
        print(f"\n  Flux stats:")
        print(f"    q_x range: [{q_x.min():.6e}, {q_x.max():.6e}]")
        print(f"    q_y range: [{q_y.min():.6e}, {q_y.max():.6e}]")
        print(f"  Omega stats:")
        print(f"    range: [{omega.min():.6e}, {omega.max():.6e}]")
        print(f"  Forces stats:")
        print(f"    fx range: [{f[:, 0].min():.6e}, {f[:, 0].max():.6e}]")
        print(f"    fy range: [{f[:, 1].min():.6e}, {f[:, 1].max():.6e}]")


def compare_restarts(file1, file2, tol=1e-10):
    """
    Compare two restart files.
    
    Args:
        file1, file2: paths to restart files
        tol: relative tolerance for comparison
    
    Returns:
        dict with comparison results
    """
    r1 = IBPMRestart(file1)
    r2 = IBPMRestart(file2)
    
    results = {
        'match': True,
        'differences': {}
    }
    
    # Check grid parameters
    for key in ['nx', 'ny', 'ngrid', 'dx', 'x0', 'y0', 'numPoints']:
        if r1.data[key] != r2.data[key]:
            results['differences'][f'grid_{key}'] = (r1.data[key], r2.data[key])
            results['match'] = False
    
    # Check timestep and time
    if r1.data['timestep'] != r2.data['timestep']:
        results['differences']['timestep'] = (r1.data['timestep'], r2.data['timestep'])
        results['match'] = False
    if abs(r1.data['time'] - r2.data['time']) > tol:
        results['differences']['time'] = (r1.data['time'], r2.data['time'])
        results['match'] = False
    
    # Check data arrays
    q_x1, q_y1 = r1.get_flux_xy()
    q_x2, q_y2 = r2.get_flux_xy()
    
    if not np.allclose(q_x1, q_x2, rtol=tol, atol=0):
        max_rel_err = np.max(np.abs((q_x1 - q_x2) / (np.abs(q_x1) + 1e-15)))
        results['differences']['q_x'] = f"max_rel_error={max_rel_err}"
        results['match'] = False
    
    if not np.allclose(q_y1, q_y2, rtol=tol, atol=0):
        max_rel_err = np.max(np.abs((q_y1 - q_y2) / (np.abs(q_y1) + 1e-15)))
        results['differences']['q_y'] = f"max_rel_error={max_rel_err}"
        results['match'] = False
    
    omega1 = r1.get_omega()
    omega2 = r2.get_omega()
    
    if not np.allclose(omega1, omega2, rtol=tol, atol=0):
        max_rel_err = np.max(np.abs((omega1 - omega2) / (np.abs(omega1) + 1e-15)))
        results['differences']['omega'] = f"max_rel_error={max_rel_err}"
        results['match'] = False
    
    f1 = r1.get_forces()
    f2 = r2.get_forces()
    
    if not np.allclose(f1, f2, rtol=tol, atol=0):
        max_rel_err = np.max(np.abs((f1 - f2) / (np.abs(f1) + 1e-15)))
        results['differences']['forces'] = f"max_rel_error={max_rel_err}"
        results['match'] = False
    
    return results


if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python validate_restart_reader.py <restart_file> [<restart_file2> ...]")
        print("\nExample: python validate_restart_reader.py cylinder_test00000.bin")
        sys.exit(1)
    
    # Read and display info for each file
    for filename in sys.argv[1:]:
        if Path(filename).exists():
            print("\n" + "=" * 70)
            r = IBPMRestart(filename)
            r.info()
        else:
            print(f"File not found: {filename}")
    
    # Compare first two files if provided
    if len(sys.argv) >= 3:
        print("\n" + "=" * 70)
        print(f"Comparing {sys.argv[1]} and {sys.argv[2]}")
        result = compare_restarts(sys.argv[1], sys.argv[2])
        print(f"Match: {result['match']}")
        if result['differences']:
            print("Differences:")
            for key, val in result['differences'].items():
                print(f"  {key}: {val}")


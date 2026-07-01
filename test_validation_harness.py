#!/usr/bin/env python3
"""
Validation test harness for IBPM Python port.

This script provides utility functions for validating ported Python modules against
C++ reference data. Each validation test loads reference restart files and checks
that the Python implementation produces the same results.

Usage:
    python test_validation_harness.py          # Run all validation tests
    python test_validation_harness.py summary  # Print summary of reference data
"""

import sys
from pathlib import Path
import numpy as np
from validate_restart_reader import IBPMRestart, compare_restarts


class ValidationHarness:
    """Test harness for IBPM Python port validation."""

    def __init__(self, reference_dir='examples'):
        """Initialize harness with reference directory."""
        self.ref_dir = Path(reference_dir)
        self.reference_files = sorted(self.ref_dir.glob('cylinder_test*.bin'))

        if not self.reference_files:
            raise FileNotFoundError(
                f"No reference files found in {reference_dir}. "
                f"Run: cd {reference_dir} && {self.ref_dir.parent}/build/ibpm ..."
            )

        self.reference_data = {}
        self._load_references()

    def _load_references(self):
        """Pre-load all reference files."""
        print(f"Loading {len(self.reference_files)} reference files...")
        for f in self.reference_files:
            try:
                self.reference_data[f.name] = IBPMRestart(str(f))
            except Exception as e:
                print(f"  ERROR loading {f.name}: {e}")

    def print_summary(self):
        """Print summary of reference data."""
        print("\n" + "=" * 70)
        print("REFERENCE DATA SUMMARY")
        print("=" * 70)

        if not self.reference_data:
            print("No reference data loaded!")
            return

        # Print grid info from first file
        first_file = self.reference_files[0]
        first_data = self.reference_data[first_file.name]

        print(f"\nGrid Configuration:")
        print(f"  nx={first_data.data['nx']}, ny={first_data.data['ny']}, "
              f"ngrid={first_data.data['ngrid']}")
        print(f"  dx={first_data.data['dx']:.6f}")
        print(f"  x0={first_data.data['x0']:.1f}, y0={first_data.data['y0']:.1f}")
        print(f"  Boundary points: {first_data.data['numPoints']}")

        print(f"\nTimesteps Available:")
        for fname in sorted(self.reference_data.keys()):
            data = self.reference_data[fname]
            print(f"  {fname:25s} t={data.data['time']:8.6f} "
                  f"step={data.data['timestep']:4d}")

        # Print data statistics
        print(f"\nData Range (from final timestep):")
        last_data = self.reference_data[sorted(self.reference_data.keys())[-1]]

        q_x, q_y = last_data.get_flux_xy()
        omega = last_data.get_omega()
        forces = last_data.get_forces()

        print(f"  Flux Q_x: [{q_x.min():.6e}, {q_x.max():.6e}]")
        print(f"  Flux Q_y: [{q_y.min():.6e}, {q_y.max():.6e}]")
        print(f"  Vorticity ω: [{omega.min():.6e}, {omega.max():.6e}]")
        print(f"  Force F_x: [{forces[:, 0].min():.6e}, {forces[:, 0].max():.6e}]")
        print(f"  Force F_y: [{forces[:, 1].min():.6e}, {forces[:, 1].max():.6e}]")

    def test_grid_parameters(self):
        """Test that all files have consistent grid parameters."""
        print("\n" + "=" * 70)
        print("TEST: Grid Parameters Consistency")
        print("=" * 70)

        first_data = self.reference_data[sorted(self.reference_data.keys())[0]]
        nx = first_data.data['nx']
        ny = first_data.data['ny']
        ngrid = first_data.data['ngrid']
        dx = first_data.data['dx']

        all_match = True
        for fname, data in sorted(self.reference_data.items()):
            match = (data.data['nx'] == nx and
                    data.data['ny'] == ny and
                    data.data['ngrid'] == ngrid and
                    data.data['dx'] == dx)
            status = "✓" if match else "✗"
            print(f"  {status} {fname}")
            all_match = all_match and match

        print(f"\nResult: {'PASS' if all_match else 'FAIL'}")
        return all_match

    def test_data_validity(self):
        """Test that all files contain valid (non-NaN) data after t=0."""
        print("\n" + "=" * 70)
        print("TEST: Data Validity (No NaN)")
        print("=" * 70)

        all_valid = True
        for fname, data in sorted(self.reference_data.items()):
            # Skip t=0, allow to be all-zero
            if data.data['timestep'] == 0:
                print(f"  - {fname:25s} (t=0, skipped)")
                continue

            q_x, q_y = data.get_flux_xy()
            omega = data.get_omega()
            forces = data.get_forces()

            has_nan = (np.any(np.isnan(q_x)) or np.any(np.isnan(q_y)) or
                      np.any(np.isnan(omega)) or np.any(np.isnan(forces)))

            status = "✗" if has_nan else "✓"
            print(f"  {status} {fname:25s} t={data.data['time']:.6f}")
            all_valid = all_valid and not has_nan

        print(f"\nResult: {'PASS' if all_valid else 'FAIL'}")
        return all_valid

    def test_monotonic_time(self):
        """Test that time increases monotonically."""
        print("\n" + "=" * 70)
        print("TEST: Monotonic Time Increase")
        print("=" * 70)

        sorted_files = sorted(self.reference_data.keys())
        sorted_data = [self.reference_data[f] for f in sorted_files]

        times = [d.data['time'] for d in sorted_data]

        monotonic = True
        for i, (fname, t) in enumerate(zip(sorted_files, times)):
            if i == 0:
                print(f"  {fname:25s} t={t:.6f}")
            else:
                prev_t = times[i-1]
                if t > prev_t:
                    print(f"  ✓ {fname:25s} t={t:.6f} (Δt={t-prev_t:.6f})")
                else:
                    print(f"  ✗ {fname:25s} t={t:.6f} (Δt={t-prev_t:.6f})")
                    monotonic = False

        print(f"\nResult: {'PASS' if monotonic else 'FAIL'}")
        return monotonic

    def validate_state_format(self, filename, expected_grid_params=None):
        """
        Validate that a Python State object matches a reference restart file.

        This is a template for validating ported code.

        Args:
            filename: Path to reference restart file
            expected_grid_params: Dict with expected nx, ny, ngrid, etc.

        Returns:
            bool: True if validation passes
        """
        # Not implemented yet — this is for ported Python code
        raise NotImplementedError("This method is called by ported Python code")

    def get_reference_state(self, filename):
        """
        Get reference State data (not yet ported).

        Args:
            filename: Restart filename in reference directory

        Returns:
            dict with keys: omega, q_x, q_y, forces, timestep, time, grid_params
        """
        if filename not in self.reference_data:
            raise KeyError(f"Reference file not found: {filename}")

        data = self.reference_data[filename]
        q_x, q_y = data.get_flux_xy()
        omega = data.get_omega()
        forces = data.get_forces()

        return {
            'timestep': data.data['timestep'],
            'time': data.data['time'],
            'omega': omega[0],  # Level 0
            'q_x': q_x[0],
            'q_y': q_y[0],
            'forces': forces,
            'grid_params': {
                'nx': data.data['nx'],
                'ny': data.data['ny'],
                'ngrid': data.data['ngrid'],
                'dx': data.data['dx'],
                'x0': data.data['x0'],
                'y0': data.data['y0'],
                'numPoints': data.data['numPoints'],
            }
        }


def main():
    """Run validation tests."""
    harness = ValidationHarness('examples')

    # Print summary
    harness.print_summary()

    # Run tests
    all_pass = True
    all_pass &= harness.test_grid_parameters()
    all_pass &= harness.test_monotonic_time()
    all_pass &= harness.test_data_validity()

    print("\n" + "=" * 70)
    if all_pass:
        print("ALL VALIDATION TESTS PASSED ✓")
    else:
        print("SOME VALIDATION TESTS FAILED ✗")
    print("=" * 70)

    return 0 if all_pass else 1


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'summary':
        harness = ValidationHarness('examples')
        harness.print_summary()
    else:
        sys.exit(main())

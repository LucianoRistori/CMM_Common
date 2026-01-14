# common — Shared CMM Utilities

This repository provides shared C++ utilities used across multiple
Coordinate Measurement Machine (CMM) analysis tools (e.g. CompareScan,
FlatnessScan, ExpansionScan, etc.).

The focus is on robust, numerically stable geometry handling and reusable
infrastructure for point-based measurements.

===============================================================================

POINTS

The Points module defines a lightweight container and I/O helpers for CMM
point data.

Each point consists of:
  - a string label
  - a vector of coordinates (coords)

The first three coordinates are interpreted as X, Y, Z in millimeters.
Additional coordinates (if present) are preserved and propagated unchanged.

Features:
  - CSV-style input with flexible column count
  - Consistent interpretation of XYZ coordinates
  - Preservation of extra per-point data
  - Shared implementation across all CMM tools

This module forms the foundation for all point-based operations in the CMM
workflow.

===============================================================================

ROTOTRANSLATION

The RotoTranslation module implements rigid 3D transformations (rotation
plus translation) in a safe, explicit, and reusable form.

A roto-translation is defined as:

  p' = R · p + t

where:
  - R is a 3×3 orthonormal rotation matrix
  - t is a 3D translation vector

This convention is used consistently throughout the code.

Design principles:

  - The rotation matrix R is the authoritative representation.
    It is always stored, applied, and inverted directly.

  - Euler angles (roll, pitch, yaw) are derived quantities.
    They are provided only for input/output and human-readable diagnostics,
    and are never used internally to apply transformations.

  - The inverse transform is exact (up to floating-point precision):

        R⁻¹ = Rᵀ
        t⁻¹ = −Rᵀ · t

Features:

  - Apply a transform to a single Point or to a vector<Point>
  - In-place or copy-based application
  - Exact inverse via matrix transpose
  - Composition of transforms
  - Euler ZYX helpers (roll–pitch–yaw), using the convention:

        R = Rz(yaw) · Ry(pitch) · Rx(roll)

  - Preservation of point labels and extra coordinates

This abstraction eliminates common errors related to sign conventions,
application order, and inverse transforms.

===============================================================================

FITROTOTRANSLATION (RIGID ALIGNMENT)

The FitRotoTranslation module computes the best rigid transformation between
two corresponding 3D point clouds.

The fit uses Horn’s closed-form quaternion solution for absolute orientation:

  B. K. P. Horn,
  "Closed-form solution of absolute orientation using unit quaternions",
  Journal of the Optical Society of America A, Vol. 4, No. 4,
  pp. 629–642 (1987).
  DOI: 10.1364/JOSAA.4.000629

Properties:

  - Least-squares optimal rigid fit
  - Closed-form (non-iterative) solution
  - Produces an orthonormal rotation matrix with determinant +1
  - No reflections
  - Numerically stable

The result is returned as a RotoTranslation object, ensuring consistent and
unambiguous use of the authoritative (R, t) representation.

Output:

  - Best-fit RotoTranslation
  - RMS residual (in the same units as input coordinates)
  - Number of points used in the fit

===============================================================================

TESTING

Two complementary tests validate the implementation:

1) Noiseless test
   - Synthetic data with an exact known transform
   - Validates correctness at machine precision (~1e−12)

2) Noisy realistic test
   - Random point cloud
   - Random rigid transform
   - Gaussian noise added independently to all coordinates
   - Fit residuals analyzed via histograms of X, Y, Z differences

These tests validate both mathematical correctness and statistical behavior
under realistic measurement noise.

===============================================================================

INTENDED USAGE

This repository is intended to be used as a shared dependency by multiple
CMM analysis tools, providing:

  - a single, trusted implementation of rigid transformations
  - consistent point handling
  - reproducible numerical behavior across projects

//==============================================================================
// File: MapWarm2Cold.h
//
// Purpose:
//   Provide a shared geometric utility to map XY coordinates measured
//   in a "warm" configuration to the corresponding predicted XY coordinates
//   in a "cold" configuration, using a 2D similarity transform.
//
//   The transform parameters (tx, ty, s, theta) are typically obtained from
//   a global fit (e.g. ExpansionScan) that aligns warm and cold measurements
//   by minimizing residuals in the XY plane.
//
//   This function implements the forward mapping:
//
//       [xc] = [tx] + s * R(theta) * [xw]
//       [yc]   [ty]
//
//   where:
//     - (xw, yw) are the warm coordinates
//     - (xc, yc) are the predicted cold coordinates
//     - s        is a uniform (isotropic) scale factor
//     - theta    is a rotation angle in radians (counter-clockwise)
//     - R(theta) is the standard 2×2 rotation matrix
//
//   The function is:
//     • Pure (no side effects)
//     • Independent of ROOT
//     • Suitable for reuse across analysis tools
//
//   Typical uses:
//     • Predict cold coordinates from warm measurements
//     • Compute per-point residuals after a global fit
//     • Share a single, consistent geometric mapping across projects
//
// Interface:
//   MapWarm2Cold(xw, yw, tx, ty, s, theta)
//
// Returns:
//   A std::pair<double,double> containing (xc, yc).
//
// Notes:
//   - Angles are in radians.
//   - All coordinates are assumed to be in the same linear units (e.g. mm).
//   - This is a forward map only; inverse mapping may be added in the future.
//
//==============================================================================

#pragma once

#include <utility>

std::pair<double,double>
MapWarm2Cold(double xw, double yw,
             double tx, double ty,
             double s, double theta);

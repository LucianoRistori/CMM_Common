//======================================================================
// File: Points.h
// Author: Luciano Ristori
// Created: original version (FlatnessScan / CompareScan common module)
// Updated: October 2025
//
// Description:
//   Defines the Point structure and the readPoints() function used by
//   FlatnessScan, CompareScan, OptimizePath, and other CMM data tools.
//
//   Each Point represents a 3D measurement (X, Y, Z) optionally labeled
//   by an identifying string (e.g. P1, A03, etc.). The label field is
//   preserved when present in the input file but left empty otherwise.
//
//   The readPoints() function reads a list of points from a text or CSV
//   file. It supports the following formats:
//
//       <label> <X> <Y> <Z>
//       <X> <Y> <Z>
//       <label>,<X>,<Y>,<Z>
//       <X>,<Y>,<Z>
//
//   Lines that cannot be parsed are skipped with a warning.
//   This version is backward compatible with previous releases.
//
// Revision history:
//   v1.2.0  (Oct 2025)  Added 'label' field to Point struct.
//                        Updated readPoints() to handle labels and CSV.
//======================================================================

#ifndef POINTS_H
#define POINTS_H

#include <string>
#include <vector>

//---------------------------------------------------------------
// Basic data structure for a measured 3D point.
//---------------------------------------------------------------
struct Point {
    std::string label;   // optional point label (empty if not present)
    double coords[3];    // X, Y, Z coordinates
};

//---------------------------------------------------------------
// Read list of points from file (CSV or space-separated).
// Returns vector of valid points. Lines that cannot be parsed
// are skipped with a warning.
//
// Parameters:
//   filename   → input file name
//   nExpected  → optional expected number of points (0 = ignore)
//
//---------------------------------------------------------------
std::vector<Point> readPoints(const std::string& filename, int nExpected = 0);

#endif // POINTS_H

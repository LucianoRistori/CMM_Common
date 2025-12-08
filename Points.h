//======================================================================
// File: Points.h
// Author: Luciano Ristori
// Created: original version (FlatnessScan / CompareScan common module)
// Updated: December 2025 (v2.0)
//
// Description:
//   Defines the Point structure and the readPoints() function used by
//   FlatnessScan, CompareScan, OptimizePath, AddDisplacedPoints, and other
//   CMM data tools.
//
//   A Point consists of:
//       label  → mandatory non-numeric identifier (string)
//       coords → a vector<double> of N coordinates, where N is specified
//                 by the caller via readPoints(filename, nExpected)
//
//   The readPoints() function enforces the following rules:
//
//     • First token MUST be a non-numeric label
//         - If the first token can be parsed as a number, the line is rejected.
//
//     • After the label, exactly nExpected numeric coordinates are required.
//         - If fewer are present → skip line with warning.
//         - If more are present  → ignore extras.
//
//     • Commas are treated as spaces (CSV compatibility).
//
//     • Empty / blank lines produce warnings and are skipped.
//
//     • Warnings use ROOT-style formatting:
//         *** Warning in readPoints(): <file> line <N>: <reason>
//         --> <offending line content>
//
//     • Trailing garbage after the coordinates is allowed and ignored.
//
//     • All warnings are printed to std::cerr.
//       The parser never aborts; it continues through the file.
//
//----------------------------------------------------------------------

#ifndef POINTS_H
#define POINTS_H

#include <string>
#include <vector>

//---------------------------------------------------------------
// Basic data structure for a measured point with arbitrary
// number of coordinates.
//---------------------------------------------------------------
struct Point {
    std::string label;           // mandatory non-numeric point label
    std::vector<double> coords;  // variable-length coordinate list
};

//---------------------------------------------------------------
// Read points from text or CSV file.
//
// Parameters:
//   filename   → name of the input file
//   nExpected  → number of coordinates expected per point
//                (must be ≥ 1)
//
// Returns:
//   vector<Point> containing all valid points.
//   Invalid lines are skipped with ROOT-style warnings.
//
//---------------------------------------------------------------
std::vector<Point> readPoints(const std::string& filename,
                              int nExpected);

#endif // POINTS_H

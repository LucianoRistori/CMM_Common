#ifndef POINTS_H
#define POINTS_H


//------------------------------------------------------------------------------
// Read a set of points from an ASCII file.
//
// Each line of the input file must contain:
//     <label> <coord1> <coord2> ... <coordN>
//
// The first field ("label") can be numeric (e.g. 1, 2, 3) or alphanumeric
// (e.g. P1, P2, P3). It is ignored during parsing.
//
// Parameters:
//   filename    - name of the input text file
//   nExpected   - number of coordinate fields to read from each line
//                 (e.g. 3 for X,Y,Z; 4 if an additional value follows).
//
// Behavior:
//   • Lines starting with '#' or empty lines are ignored.
//   • Commas are treated as spaces (CSV compatible).
//   • Lines with fewer than nExpected numeric values are skipped with a warning.
//   • Lines with extra numbers are accepted but only the first nExpected
//     coordinates are stored.
//   • The first field (label) is ignored; each Point receives an internal
//     sequential ID starting from 0.
//   • Returns a vector of valid Point structures.
//
// Example input:
//     P1 12.345 67.890 0.123
//     P2 12.346 67.891 0.122
//
// Example call:
//     auto points = readPoints("data.txt", 3);
//
// Example output (2 valid points read):
//     Warning: line 3 has only 2 numeric fields (expected at least 3). Skipped.
//     Read 2 valid points from data.txt
//------------------------------------------------------------------------------
#include <string>
#include <vector>

//------------------------------------------------------------------------------
// Structure representing a single measurement point
//------------------------------------------------------------------------------
struct Point {
    int id = 0;                  // Sequential internal ID assigned by readPoints()
    double coords[10] = {0.0};   // Coordinate array (size should be >= nExpected)
};

std::vector<Point> readPoints(const std::string &filename, int nExpected);

#endif // POINTS_H

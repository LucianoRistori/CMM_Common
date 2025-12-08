#include "Points.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>

using namespace std;

//--------------------------------------------------------------------------
// Helper: return true if token is numeric (integer or floating-point).
//--------------------------------------------------------------------------
static bool isNumericToken(const string& s) {
    if (s.empty()) return false;

    // Allow leading + or -
    size_t i = 0;
    if (s[i] == '+' || s[i] == '-') i++;

    bool hasDigit = false;
    bool hasDot   = false;

    for (; i < s.size(); ++i) {
        char c = s[i];
        if (isdigit(c)) {
            hasDigit = true;
            continue;
        }
        if (c == '.') {
            if (hasDot) return false; // second dot → invalid number
            hasDot = true;
            continue;
        }
        // Any other character makes it non-numeric
        return false;
    }

    return hasDigit; // e.g. "." alone is NOT numeric
}

//--------------------------------------------------------------------------
// readPoints() – Version 3.1
//
// Requirements:
//   • First token must be a NON-NUMERIC label (e.g. P1, A03, G12)
//   • If first meaningful data line has numeric first token → FATAL ERROR
//   • Commas allowed (converted to spaces)
//   • Blank lines produce warnings
//   • Exactly nExpected coordinates must be provided per line
//     (additional tokens ignored; missing tokens = skip line with warning)
//
// Behavior:
//   • Stops reading and aborts at the FIRST malformed label line.
//   • Does NOT return partial data in such a case.
//--------------------------------------------------------------------------
vector<Point> readPoints(const string& filename, int nExpected) {

    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "*** ERROR in readPoints(): cannot open file \""
             << filename << "\"\n";
        exit(1);
    }

    vector<Point> points;
    string line;
    int lineNum = 0;
    bool firstPointSeen = false;

    while (getline(infile, line)) {
        ++lineNum;

        // Trim whitespace-only lines
        bool blank = true;
        for (char c : line) {
            if (!isspace((unsigned char)c)) { blank = false; break; }
        }
        if (blank) {
            cerr << "*** Warning: " << filename << " line " << lineNum
                 << ": empty or blank line\n";
            continue;
        }

        // Replace commas → spaces
        for (char &c : line) {
            if (c == ',') c = ' ';
        }

        // Tokenize
        istringstream iss(line);
        vector<string> tokens;
        string tok;

        while (iss >> tok)
            tokens.push_back(tok);

        if (tokens.empty()) {
            cerr << "*** Warning: " << filename << " line " << lineNum
                 << ": no usable tokens\n";
            continue;
        }

        //------------------------------------------------------------------
        // FIRST meaningful point line: enforce NON-numeric label rule.
        //------------------------------------------------------------------
        if (!firstPointSeen) {
            const string& labelToken = tokens[0];

            if (isNumericToken(labelToken)) {
                cerr << "\n*** FATAL ERROR in readPoints(\"" << filename << "\"):\n"
                     << "    Line " << lineNum << ": first token is numeric (\""
                     << labelToken << "\"), but a non-numeric label is required.\n"
                     << "    Every point line must begin with a string label.\n"
                     << "    Example:  P01  10.0  20.0  30.0\n"
                     << "              A3   12.5  18.0  29.7\n"
                     << "    Aborting.\n\n";
                exit(1);
            }

            firstPointSeen = true;
        }

        //------------------------------------------------------------------
        // Extract label and coordinates
        //------------------------------------------------------------------
        Point p;
        p.label = tokens[0];

        if (nExpected <= 0) {
            cerr << "*** ERROR in readPoints(): nExpected must be > 0\n";
            exit(1);
        }

        p.coords.clear();
        p.coords.reserve(nExpected);

        // Need at least 1 label + nExpected coordinates
        if ((int)tokens.size() < (1 + nExpected)) {
            cerr << "*** Warning: " << filename << " line " << lineNum
                 << ": expected " << nExpected << " coordinates after label \""
                 << p.label << "\", but found only " << (tokens.size() - 1)
                 << ". Skipping line.\n";
            continue;
        }

        // Parse exactly nExpected coordinates
        bool bad = false;
        for (int i = 0; i < nExpected; ++i) {
            double v = 0.0;
            const string& s = tokens[1 + i];
            try {
                size_t idx = 0;
                v = stod(s, &idx);
                if (idx != s.size()) {
                    bad = true;
                    break;
                }
            } catch (...) {
                bad = true;
                break;
            }
            p.coords.push_back(v);
        }

        if (bad) {
            cerr << "*** Warning: " << filename << " line " << lineNum
                 << ": invalid numeric value in coordinates after label \""
                 << p.label << "\". Skipping line.\n";
            continue;
        }

        //------------------------------------------------------------------
        // Append the parsed point
        //------------------------------------------------------------------
        points.push_back(p);
    }

    return points;
}

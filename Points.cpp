#include "Points.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

vector<Point> readPoints(const string &filename, int nExpected) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        exit(1);
    }

    vector<Point> points;
    string line;
    int lineNum = 0;

    while (getline(infile, line)) {
        ++lineNum;
        if (line.empty()) continue;

        // Replace commas with spaces for CSV files
        for (char &c : line)
            if (c == ',') c = ' ';

        istringstream iss(line);
        Point p;
        p.label = "";

        // Try to read label + 3 numbers
        if (!(iss >> p.label >> p.coords[0] >> p.coords[1] >> p.coords[2])) {
            // Retry without label (3 numeric columns)
            iss.clear();
            iss.str(line);
            if (!(iss >> p.coords[0] >> p.coords[1] >> p.coords[2])) {
                cerr << "Warning: skipping invalid line " << lineNum << endl;
                continue;
            }
        }

        points.push_back(p);
    }

    if (nExpected > 0 && (int)points.size() != nExpected)
        cerr << "Warning: expected " << nExpected
             << " points but read " << points.size() << endl;

    return points;
}

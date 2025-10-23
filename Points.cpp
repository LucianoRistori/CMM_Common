#include "Points.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cctype>

std::vector<Point> readPoints(const std::string &filename, int nExpected) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(1);
    }

    std::vector<Point> points;
    std::string line;
    int lineNum = 0;
    int index = 0;  // sequential internal ID

    while (std::getline(infile, line)) {
        ++lineNum;

        // Replace commas with spaces (allow CSV or space-delimited)
        for (char &c : line) {
            if (c == ',') c = ' ';
        }

        // Skip empty or comment lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string label;
        iss >> label;  // skip first field (label)

        // Read numeric coordinates
        std::vector<double> coords;
        double val;
        while (iss >> val) coords.push_back(val);

        if ((int)coords.size() < nExpected) {
            std::cout << "Warning: line " << lineNum
                      << " has only " << coords.size()
                      << " numeric fields (expected at least "
                      << nExpected << "). Skipped.\n";
            continue;
        }

        // Create Point and fill only the first nExpected coordinates
        Point p;
        p.id = index++;
        for (int i = 0; i < nExpected; ++i) {
            p.coords[i] = coords[i];
        }

        points.push_back(p);
    }

    infile.close();

    if (points.empty()) {
        std::cerr << "Error: no valid points read from file " << filename << std::endl;
    } else {
        std::cout << "Read " << points.size()
                  << " valid points from " << filename << std::endl;
    }

    return points;
}

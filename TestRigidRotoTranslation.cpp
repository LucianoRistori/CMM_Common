//==============================================================================
// File: TestRigidRotoTranslation.cpp
//==============================================================================
//
// Build example (from common directory):
// clang++ -std=c++17 -O2 -Wall -Wextra TestRigidRotoTranslation.cpp Points.cpp -I. -o TestRigidRotoTranslation
//
// (Adjust include path / Points.cpp path as needed for your layout.)
//==============================================================================

#include <iostream>
#include <iomanip>
#include <cmath>
#include "RigidRotoTranslation.h"

static void printPts(const std::vector<Point>& v, const char* title){
    std::cout << title << "\n";
    for(const auto& p : v){
        std::cout << "  " << std::setw(8) << p.label << " : "
                  << "(" << p.coords[0] << ", " << p.coords[1] << ", " << p.coords[2] << ")";
        if(p.coords.size() > 3){
            std::cout << "  [extra coords preserved: " << (p.coords.size()-3) << "]";
        }
        std::cout << "\n";
    }
}

int main(){
    // Create a tiny point set compatible with Points.h::Point
    std::vector<Point> pts;
    {
        Point a; a.label="A"; a.coords={1.0, 0.0, 0.0};
        Point b; b.label="B"; b.coords={0.0, 1.0, 0.0};
        Point c; c.label="C"; c.coords={0.0, 0.0, 1.0, 123.0}; // extra coord to prove we preserve it
        pts.push_back(a); pts.push_back(b); pts.push_back(c);
    }

    printPts(pts, "Original points:");

    const double deg = M_PI / 180.0;

    // Test 1: pure translation
    auto t1 = RigidRotoTranslation(pts, 10.0, 20.0, 30.0, 0.0, 0.0, 0.0);
    printPts(t1, "\nTest 1: translation (10,20,30), no rotation:");

    // Test 2: 90 deg yaw about Z, no translation
    auto t2 = RigidRotoTranslation(pts, 0.0, 0.0, 0.0, 0.0, 0.0, 90.0*deg);
    printPts(t2, "\nTest 2: 90 deg rotation about Z (yaw), no translation:");
    std::cout << "Expected (approximately):\n"
              << "  A(1,0,0) -> (0,1,0)\n"
              << "  B(0,1,0) -> (-1,0,0)\n";

    // Test 3: rotation + translation (in-place)
    auto t3 = pts;
    RigidRotoTranslationInPlace(t3, 1.0, 2.0, 3.0, 0.0, 0.0, 90.0*deg);
    printPts(t3, "\nTest 3: 90 deg Z rotation + translation (1,2,3), in-place:");

    return 0;
}
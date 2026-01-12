//==============================================================================
// File: RigidRotoTranslation.h
// Purpose:
//   Apply a rigid 6-DOF transform (rotation + translation) to points.
//
// Integration with common/Points.h:
//   This header provides overloads that operate directly on the Point type
//   defined in Points.h. It rotates/translates the first 3 coordinates
//   (X,Y,Z) and preserves:
//     • label (if present in Point)
//     • any extra coordinates beyond the first 3 (if present)
//
// Convention:
//   p' = Rz(yaw) * Ry(pitch) * Rx(roll) * p + T
//
// Angles are in radians.
//==============================================================================

#ifndef RIGID_ROTOTRANSLATION_H
#define RIGID_ROTOTRANSLATION_H

#include <vector>
#include <cmath>
#include <stdexcept>

#include "Points.h"

// --- Helper: apply Rz*Ry*Rx to a single (x,y,z) and add translation ----------
inline void ApplyRigidRotoTranslationXYZ(
    double x, double y, double z,
    double tx, double ty, double tz,
    double roll, double pitch, double yaw,
    double& xo, double& yo, double& zo
){
    const double cr = std::cos(roll),  sr = std::sin(roll);
    const double cp = std::cos(pitch), sp = std::sin(pitch);
    const double cy = std::cos(yaw),   sy = std::sin(yaw);

    // Rotation matrix: R = Rz(yaw) * Ry(pitch) * Rx(roll)
    const double r00 =  cy * cp;
    const double r01 =  cy * sp * sr - sy * cr;
    const double r02 =  cy * sp * cr + sy * sr;

    const double r10 =  sy * cp;
    const double r11 =  sy * sp * sr + cy * cr;
    const double r12 =  sy * sp * cr - cy * sr;

    const double r20 = -sp;
    const double r21 =  cp * sr;
    const double r22 =  cp * cr;

    xo = r00*x + r01*y + r02*z + tx;
    yo = r10*x + r11*y + r12*z + ty;
    zo = r20*x + r21*y + r22*z + tz;
}

// --- Overload: transform a vector<Point> (from Points.h) ---------------------
inline std::vector<Point> RigidRotoTranslation(
    const std::vector<Point>& pts,
    double tx, double ty, double tz,
    double roll, double pitch, double yaw
){
    std::vector<Point> out;
    out.reserve(pts.size());

    for(const auto& p : pts){
        if(p.coords.size() < 3){
            throw std::runtime_error("RigidRotoTranslation: Point has < 3 coordinates.");
        }

        Point q = p; // preserves label + any extra coords
        double xo=0, yo=0, zo=0;

        ApplyRigidRotoTranslationXYZ(
            p.coords[0], p.coords[1], p.coords[2],
            tx, ty, tz, roll, pitch, yaw,
            xo, yo, zo
        );

        q.coords[0] = xo;
        q.coords[1] = yo;
        q.coords[2] = zo;

        out.push_back(q);
    }

    return out;
}

// --- In-place variant -------------------------------------------------------
inline void RigidRotoTranslationInPlace(
    std::vector<Point>& pts,
    double tx, double ty, double tz,
    double roll, double pitch, double yaw
){
    for(auto& p : pts){
        if(p.coords.size() < 3){
            throw std::runtime_error("RigidRotoTranslationInPlace: Point has < 3 coordinates.");
        }

        double xo=0, yo=0, zo=0;
        ApplyRigidRotoTranslationXYZ(
            p.coords[0], p.coords[1], p.coords[2],
            tx, ty, tz, roll, pitch, yaw,
            xo, yo, zo
        );

        p.coords[0] = xo;
        p.coords[1] = yo;
        p.coords[2] = zo;
    }
}

#endif // RIGID_ROTOTRANSLATION_H

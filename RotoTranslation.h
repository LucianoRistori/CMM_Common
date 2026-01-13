//==============================================================================
// File: RotoTranslation.h
//
// Purpose:
//   Represent and apply a rigid 3D transform:
//
//     p' = R * p + t
//
//   where R is 3x3 orthonormal rotation and t is translation vector.
//
// Features:
//   • Apply to Point / vector<Point> (Points.h)
//   • Inverse (uses R^T, t' = -R^T t)
//   • Composition (A∘B): p' = A(B(p))
//   • Euler helpers (ZYX): R = Rz(yaw)*Ry(pitch)*Rx(roll)
//
// Notes:
//   R is authoritative. Euler angles are for I/O only.
//==============================================================================

#ifndef ROTOTRANSLATION_H
#define ROTOTRANSLATION_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "Points.h"

struct RotoTranslation {
    // Authoritative representation
    double R[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
    double t[3]    = { 0,0,0 };

    // --- Apply to raw xyz ----------------------------------------------------
    inline void apply(double x, double y, double z,
                      double& xo, double& yo, double& zo) const
    {
        xo = R[0][0]*x + R[0][1]*y + R[0][2]*z + t[0];
        yo = R[1][0]*x + R[1][1]*y + R[1][2]*z + t[1];
        zo = R[2][0]*x + R[2][1]*y + R[2][2]*z + t[2];
    }

    // --- Apply to a Point (preserves label + extra coords) -------------------
    inline Point apply(const Point& p) const
    {
        if(p.coords.size() < 3)
            throw std::runtime_error("RotoTranslation::apply(Point): Point has < 3 coordinates.");

        Point q = p;
        double xo, yo, zo;
        apply(p.coords[0], p.coords[1], p.coords[2], xo, yo, zo);
        q.coords[0] = xo; q.coords[1] = yo; q.coords[2] = zo;
        return q;
    }

    // --- Apply to vector<Point> ---------------------------------------------
    inline std::vector<Point> apply(const std::vector<Point>& pts) const
    {
        std::vector<Point> out;
        out.reserve(pts.size());
        for(const auto& p : pts) out.push_back(apply(p));
        return out;
    }

    // --- Inverse transform: p = R^T*(p' - t) --------------------------------
    inline RotoTranslation inverse() const
    {
        RotoTranslation inv;

        // inv.R = R^T
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                inv.R[i][j] = R[j][i];

        // inv.t = -R^T * t
        inv.t[0] = -(inv.R[0][0]*t[0] + inv.R[0][1]*t[1] + inv.R[0][2]*t[2]);
        inv.t[1] = -(inv.R[1][0]*t[0] + inv.R[1][1]*t[1] + inv.R[1][2]*t[2]);
        inv.t[2] = -(inv.R[2][0]*t[0] + inv.R[2][1]*t[1] + inv.R[2][2]*t[2]);

        return inv;
    }

    // --- Composition: A∘B means first B then A ------------------------------
    // If C = A.compose(B), then C(p) = A(B(p)).
    inline RotoTranslation compose(const RotoTranslation& B) const
    {
        const RotoTranslation& A = *this;
        RotoTranslation C;

        // C.R = A.R * B.R
        for(int i=0;i<3;++i){
            for(int j=0;j<3;++j){
                C.R[i][j] = 0.0;
                for(int k=0;k<3;++k) C.R[i][j] += A.R[i][k] * B.R[k][j];
            }
        }

        // C.t = A.R * B.t + A.t
        C.t[0] = A.R[0][0]*B.t[0] + A.R[0][1]*B.t[1] + A.R[0][2]*B.t[2] + A.t[0];
        C.t[1] = A.R[1][0]*B.t[0] + A.R[1][1]*B.t[1] + A.R[1][2]*B.t[2] + A.t[1];
        C.t[2] = A.R[2][0]*B.t[0] + A.R[2][1]*B.t[1] + A.R[2][2]*B.t[2] + A.t[2];

        return C;
    }

    // --- Euler helpers (ZYX intrinsic): R = Rz(yaw)*Ry(pitch)*Rx(roll) -------
    static inline RotoTranslation fromEulerZYX(double tx, double ty, double tz,
                                               double roll, double pitch, double yaw)
    {
        RotoTranslation T;
        T.t[0]=tx; T.t[1]=ty; T.t[2]=tz;

        const double cr = std::cos(roll),  sr = std::sin(roll);
        const double cp = std::cos(pitch), sp = std::sin(pitch);
        const double cy = std::cos(yaw),   sy = std::sin(yaw);

        T.R[0][0] =  cy*cp;
        T.R[0][1] =  cy*sp*sr - sy*cr;
        T.R[0][2] =  cy*sp*cr + sy*sr;

        T.R[1][0] =  sy*cp;
        T.R[1][1] =  sy*sp*sr + cy*cr;
        T.R[1][2] =  sy*sp*cr - cy*sr;

        T.R[2][0] = -sp;
        T.R[2][1] =  cp*sr;
        T.R[2][2] =  cp*cr;

        return T;
    }

    // Derive Euler angles from R (for printing / I-O)
    inline void toEulerZYX(double& roll, double& pitch, double& yaw) const
    {
        const double sp = -R[2][0];
        const double spc = std::max(-1.0, std::min(1.0, sp));
        pitch = std::asin(spc);

        const double cp = std::cos(pitch);

        if(std::fabs(cp) < 1e-12){
            roll = 0.0;
            yaw  = std::atan2(-R[0][1], R[1][1]);
        } else {
            roll = std::atan2(R[2][1], R[2][2]);
            yaw  = std::atan2(R[1][0], R[0][0]);
        }
    }
};

#endif // ROTOTRANSLATION_H

//==============================================================================
// File: TestFitRotoTranslation.cpp
//
// Purpose:
//   End-to-end test for:
//     • RotoTranslation (apply / inverse / compose / Euler helpers)
//     • FitRotoTranslation (Horn quaternion fit)
//
// Strategy:
//   1) Build a small non-degenerate point cloud A.
//   2) Build a "truth" transform Ttruth from Euler+T.
//   3) Generate B = Ttruth.apply(A).
//   4) Fit Tfit from A -> B.
//   5) Validate numerically using ONLY authoritative (R,t):
//        • max |Tfit(A) - B|
//        • max |Tfit^{-1}(B) - A|
//        • (optional) Tfit^{-1} ∘ Tfit ~ Identity
//
// Build (from common/):
//   clang++ -std=c++17 -O2 -Wall -Wextra TestFitRotoTranslation.cpp Points.cpp -I. -o TestFitRotoTranslation
//==============================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

#include "Points.h"
#include "RotoTranslation.h"
#include "FitRotoTranslation.h"

static double pi(){ return std::acos(-1.0); }
static double rad2deg(double r){ return r * 180.0 / pi(); }

static double wrapToPi(double a){
    const double p = pi();
    while(a >  p) a -= 2.0*p;
    while(a < -p) a += 2.0*p;
    return a;
}

static double maxAbsCoordDiff3(const std::vector<Point>& A, const std::vector<Point>& B){
    const std::size_t n = std::min(A.size(), B.size());
    double m = 0.0;
    for(std::size_t i=0;i<n;++i){
        for(int k=0;k<3;++k){
            m = std::max(m, std::fabs(A[i].coords[k] - B[i].coords[k]));
        }
    }
    return m;
}

static std::vector<Point> applyRT(const RotoTranslation& T, const std::vector<Point>& pts){
    return T.apply(pts);
}

static void printTruthFit(const RotoTranslation& Ttruth, const RotoTranslation& Tfit){
    double rT,pT,yT, rF,pF,yF;
    Ttruth.toEulerZYX(rT,pT,yT);
    Tfit.toEulerZYX  (rF,pF,yF);

    std::cout << "Truth:\n";
    std::cout << "  T = (" << Ttruth.t[0] << ", " << Ttruth.t[1] << ", " << Ttruth.t[2] << ")\n";
    std::cout << "  roll,pitch,yaw [deg] = ("
              << rad2deg(rT) << ", " << rad2deg(pT) << ", " << rad2deg(yT) << ")\n\n";

    std::cout << "Fit:\n";
    std::cout << "  T = (" << Tfit.t[0] << ", " << Tfit.t[1] << ", " << Tfit.t[2] << ")\n";
    std::cout << "  roll,pitch,yaw [deg] = ("
              << rad2deg(rF) << ", " << rad2deg(pF) << ", " << rad2deg(yF) << ")\n";
}

int main(){
    std::cout << std::fixed << std::setprecision(15);

    // -------------------------------------------------------------------------
    // 1) Build a small non-degenerate point cloud A
    // -------------------------------------------------------------------------
    std::vector<Point> A;
    auto add = [&](const char* lbl, double x, double y, double z, double extra=0.0, bool addExtra=false){
        Point p; p.label = lbl; p.coords = {x,y,z};
        if(addExtra) p.coords.push_back(extra); // ensure we preserve extra coords
        A.push_back(p);
    };

    add("P0",  0.0,  0.0,  0.0);
    add("P1", 10.0,  0.0,  0.0);
    add("P2",  0.0, 20.0,  0.0);
    add("P3",  0.0,  0.0, 30.0);
    add("P4",  5.0,  7.0, 11.0);
    add("P5", -3.0,  4.0,  9.0, 123.0, true); // add extra coord on one point

    // -------------------------------------------------------------------------
    // 2) Truth transform parameters and construction
    // -------------------------------------------------------------------------
    const double deg = pi()/180.0;

    const double tx_true =  12.3;
    const double ty_true =  -4.5;
    const double tz_true =   7.8;

    const double roll_true  =  2.0 * deg;   // X
    const double pitch_true = -3.0 * deg;   // Y
    const double yaw_true   = 15.0 * deg;   // Z

    const RotoTranslation Ttruth = RotoTranslation::fromEulerZYX(
        tx_true, ty_true, tz_true,
        roll_true, pitch_true, yaw_true
    );

    // -------------------------------------------------------------------------
    // 3) Generate B = Ttruth(A)
    // -------------------------------------------------------------------------
    const std::vector<Point> B = applyRT(Ttruth, A);

    // -------------------------------------------------------------------------
    // 4) Fit A -> B
    // -------------------------------------------------------------------------
    const FitRTResult fit = FitRotoTranslation(A, B);
    const RotoTranslation& Tfit = fit.T;

    printTruthFit(Ttruth, Tfit);

    // -------------------------------------------------------------------------
    // 5) Validate using authoritative (R,t) only
    // -------------------------------------------------------------------------
    const std::vector<Point> Bfit = applyRT(Tfit, A);
    const double maxDiff_B = maxAbsCoordDiff3(B, Bfit);

    const RotoTranslation Tinv = Tfit.inverse();
    const std::vector<Point> Aback = applyRT(Tinv, B);
    const double maxDiff_A = maxAbsCoordDiff3(A, Aback);

    std::cout << "  RMS residual = " << fit.rms << "\n\n";

    // Compare derived Euler for info only (do NOT use for validation)
    double rF,pF,yF;
    Tfit.toEulerZYX(rF,pF,yF);

    const double dtx = Tfit.t[0] - tx_true;
    const double dty = Tfit.t[1] - ty_true;
    const double dtz = Tfit.t[2] - tz_true;

    const double droll  = wrapToPi(rF - roll_true);
    const double dpitch = wrapToPi(pF - pitch_true);
    const double dyaw   = wrapToPi(yF - yaw_true);

    std::cout << "Info-only differences (Euler extraction can be touchy):\n";
    std::cout << "  dT = (" << dtx << ", " << dty << ", " << dtz << ")\n";
    std::cout << "  droll,dpitch,dyaw [deg] = ("
              << rad2deg(droll) << ", " << rad2deg(dpitch) << ", " << rad2deg(dyaw) << ")\n\n";

    std::cout << "Sanity checks (authoritative R,t):\n";
    std::cout << "  max |Tfit(A) - B|    = " << maxDiff_B << "\n";
    std::cout << "  max |Tfit^{-1}(B)-A| = " << maxDiff_A << "\n";

    // Optional: check that inverse(compose) is near identity
    const RotoTranslation I = Tinv.compose(Tfit);
    double Imax = 0.0;
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            const double target = (i==j ? 1.0 : 0.0);
            Imax = std::max(Imax, std::fabs(I.R[i][j] - target));
        }
        Imax = std::max(Imax, std::fabs(I.t[i] - 0.0));
    }
    std::cout << "  max |(Tinv∘Tfit) - I| = " << Imax << "\n\n";

    // Pass/fail thresholds (tight; should pass on synthetic noiseless data)
    const double tol = 1e-10;
    const bool ok = (fit.rms < tol) && (maxDiff_B < tol) && (maxDiff_A < tol) && (Imax < tol);

    std::cout << "Result: " << (ok ? "PASS" : "CHECK") << "\n";
    return ok ? 0 : 1;
}

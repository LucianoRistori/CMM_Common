//==============================================================================
// File: TestFitRotoTranslationNoise.cpp
//
// Purpose:
//   Realistic test of FitRotoTranslation with Gaussian noise:
//
//     1) Generate random RotoTranslation (truth)
//     2) Generate random point cloud A
//     3) Bclean = Ttruth(A)
//     4) Add Gaussian noise to Bclean -> Bnoisy
//     5) Fit Tfit from A -> Bnoisy
//     6) C = Tfit(A)
//     7) Histogram residuals (C - Bnoisy) for X, Y, Z separately
//
// Build (from common/):
//   clang++ -std=c++17 -O2 -Wall -Wextra TestFitRotoTranslationNoise.cpp Points.cpp \
//     -I. `root-config --cflags --libs` -o TestFitRotoTranslationNoise
//
// Run:
//   ./TestFitRotoTranslationNoise [Npoints=1000] [sigma=0.010] [seed=12345]
//==============================================================================

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <algorithm>

#include "Points.h"
#include "RotoTranslation.h"
#include "FitRotoTranslation.h"

// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

static double pi() { return std::acos(-1.0); }

//------------------------------------------------------------------------------
// Add Gaussian noise to XYZ only (preserve labels & extra coords)
//------------------------------------------------------------------------------
static std::vector<Point> AddGaussianNoiseXYZ(const std::vector<Point>& pts,
                                              std::mt19937_64& rng,
                                              double sigma)
{
    std::normal_distribution<double> gaus(0.0, sigma);

    std::vector<Point> out;
    out.reserve(pts.size());

    for (const auto& p : pts) {
        if (p.coords.size() < 3)
            throw std::runtime_error("Point has < 3 coordinates");

        Point q = p;
        q.coords[0] += gaus(rng);
        q.coords[1] += gaus(rng);
        q.coords[2] += gaus(rng);
        out.push_back(q);
    }
    return out;
}

//==============================================================================
// MAIN
//==============================================================================
int main(int argc, char** argv)
{
    int    N     = (argc > 1) ? std::stoi(argv[1]) : 1000;
    double sigma = (argc > 2) ? std::stod(argv[2]) : 0.010;
    uint64_t seed = (argc > 3) ? std::stoull(argv[3]) : 12345ULL;

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "Npoints = " << N << "\n";
    std::cout << "sigma   = " << sigma << "\n";
    std::cout << "seed    = " << seed << "\n\n";

    std::mt19937_64 rng(seed);

    //--------------------------------------------------------------------------
    // 1) Random truth transform
    //--------------------------------------------------------------------------
    std::uniform_real_distribution<double> uT(-50.0, 50.0);
    std::uniform_real_distribution<double> uA(-20.0*pi()/180.0,
                                                20.0*pi()/180.0);

    const double tx = uT(rng);
    const double ty = uT(rng);
    const double tz = uT(rng);

    const double roll  = uA(rng);
    const double pitch = uA(rng);
    const double yaw   = uA(rng);

    const RotoTranslation Ttruth =
        RotoTranslation::fromEulerZYX(tx, ty, tz, roll, pitch, yaw);

    //--------------------------------------------------------------------------
    // 2) Random point cloud A
    //--------------------------------------------------------------------------
    std::uniform_real_distribution<double> uX(-500.0, 500.0);

    std::vector<Point> A;
    A.reserve(N);

    for (int i = 0; i < N; ++i) {
        Point p;
        p.label = "P" + std::to_string(i);
        p.coords = { uX(rng), uX(rng), uX(rng) };
        A.push_back(p);
    }

    //--------------------------------------------------------------------------
    // 3) Bclean = Ttruth(A)
    //--------------------------------------------------------------------------
    const std::vector<Point> Bclean = Ttruth.apply(A);

    //--------------------------------------------------------------------------
    // 4) Add Gaussian noise -> Bnoisy
    //--------------------------------------------------------------------------
    const std::vector<Point> Bnoisy = AddGaussianNoiseXYZ(Bclean, rng, sigma);

    //--------------------------------------------------------------------------
    // 5) Fit A -> Bnoisy
    //--------------------------------------------------------------------------
    const FitRTResult fit = FitRotoTranslation(A, Bnoisy);
    const RotoTranslation& Tfit = fit.T;

    std::cout << "Fit RMS residual = " << fit.rms << "\n\n";

    //--------------------------------------------------------------------------
    // 6) Apply fit to A -> C
    //--------------------------------------------------------------------------
    const std::vector<Point> C = Tfit.apply(A);

    //--------------------------------------------------------------------------
    // 7) Histogram residuals: C - Bnoisy
    //--------------------------------------------------------------------------
    const double hRange = std::max(5.0 * sigma, 1e-6);
    const int nbins = 120;

    TH1D* hdx = new TH1D("residual_dx", "Residual dx = Cx - Bx;dx;Entries",
                         nbins, -hRange, +hRange);
    TH1D* hdy = new TH1D("residual_dy", "Residual dy = Cy - By;dy;Entries",
                         nbins, -hRange, +hRange);
    TH1D* hdz = new TH1D("residual_dz", "Residual dz = Cz - Bz;dz;Entries",
                         nbins, -hRange, +hRange);

    for (int i = 0; i < N; ++i) {
        hdx->Fill(C[i].coords[0] - Bnoisy[i].coords[0]);
        hdy->Fill(C[i].coords[1] - Bnoisy[i].coords[1]);
        hdz->Fill(C[i].coords[2] - Bnoisy[i].coords[2]);
    }

    //--------------------------------------------------------------------------
    // Save outputs
    //--------------------------------------------------------------------------
    TFile fout("TestFitRotoTranslationNoise.root", "RECREATE");
    hdx->Write();
    hdy->Write();
    hdz->Write();
    fout.Close();

    TCanvas c1("c1","dx",800,600); hdx->Draw(); c1.SaveAs("residual_dx.png");
    TCanvas c2("c2","dy",800,600); hdy->Draw(); c2.SaveAs("residual_dy.png");
    TCanvas c3("c3","dz",800,600); hdz->Draw(); c3.SaveAs("residual_dz.png");

    std::cout << "Histogram stats:\n";
    std::cout << "  dx: mean=" << hdx->GetMean()
              << "  rms=" << hdx->GetRMS() << "\n";
    std::cout << "  dy: mean=" << hdy->GetMean()
              << "  rms=" << hdy->GetRMS() << "\n";
    std::cout << "  dz: mean=" << hdz->GetMean()
              << "  rms=" << hdz->GetRMS() << "\n";

    return 0;
}

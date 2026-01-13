#ifndef FIT_ROTOTRANSLATION_H
#define FIT_ROTOTRANSLATION_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <cstddef>

#include "Points.h"
#include "RotoTranslation.h"

// Reference:
//   B. K. P. Horn,
//   "Closed-form solution of absolute orientation using unit quaternions,"
//   Journal of the Optical Society of America A, Vol. 4, No. 4,
//   pp. 629–642 (1987).
//   DOI: 10.1364/JOSAA.4.000629

struct FitRTResult {
    RotoTranslation T;      // authoritative transform: p' = R p + t
    double rms = 0;         // RMS residual (same units as inputs)
    std::size_t n = 0;      // points used
};

inline double Norm4(const double q[4]) {
    return std::sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
}

inline void QuaternionToR(const double q[4], double R[3][3]) {
    const double w=q[0], x=q[1], y=q[2], z=q[3];
    const double ww=w*w, xx=x*x, yy=y*y, zz=z*z;

    R[0][0] = ww + xx - yy - zz;
    R[0][1] = 2*(x*y - w*z);
    R[0][2] = 2*(x*z + w*y);

    R[1][0] = 2*(x*y + w*z);
    R[1][1] = ww - xx + yy - zz;
    R[1][2] = 2*(y*z - w*x);

    R[2][0] = 2*(x*z - w*y);
    R[2][1] = 2*(y*z + w*x);
    R[2][2] = ww - xx - yy + zz;
}

inline void LargestEigenvector4x4_Symmetric_Jacobi(const double N_in[4][4], double q_out[4])
{
    double A[4][4];
    double V[4][4] = { {1,0,0,0},
                       {0,1,0,0},
                       {0,0,1,0},
                       {0,0,0,1} };
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) A[i][j] = N_in[i][j];

    auto absd = [](double x){ return x < 0 ? -x : x; };
    const int maxSweeps = 50;
    const double eps = 1e-15;

    for(int sweep=0; sweep<maxSweeps; ++sweep){
        int p=0, r=1;
        double maxOff = absd(A[0][1]);
        for(int i=0;i<4;++i){
            for(int j=i+1;j<4;++j){
                const double v = absd(A[i][j]);
                if(v > maxOff){ maxOff=v; p=i; r=j; }
            }
        }
        if(maxOff < eps) break;

        const double app=A[p][p], arr=A[r][r], apr=A[p][r];
        const double tau = (arr - app) / (2.0 * apr);
        const double t = (tau >= 0.0)
            ? 1.0 / (tau + std::sqrt(1.0 + tau*tau))
            : -1.0 / (-tau + std::sqrt(1.0 + tau*tau));
        const double c = 1.0 / std::sqrt(1.0 + t*t);
        const double s = t * c;

        A[p][p] = app - t*apr;
        A[r][r] = arr + t*apr;
        A[p][r] = A[r][p] = 0.0;

        for(int k=0;k<4;++k){
            if(k==p || k==r) continue;
            const double akp=A[k][p], akr=A[k][r];
            A[k][p] = A[p][k] = c*akp - s*akr;
            A[k][r] = A[r][k] = c*akr + s*akp;
        }

        for(int k=0;k<4;++k){
            const double vkp=V[k][p], vkr=V[k][r];
            V[k][p] = c*vkp - s*vkr;
            V[k][r] = c*vkr + s*vkp;
        }
    }

    int imax=0;
    double lmax=A[0][0];
    for(int i=1;i<4;++i){
        if(A[i][i] > lmax){ lmax=A[i][i]; imax=i; }
    }

    q_out[0]=V[0][imax];
    q_out[1]=V[1][imax];
    q_out[2]=V[2][imax];
    q_out[3]=V[3][imax];

    const double n = Norm4(q_out);
    if(n <= 0) throw std::runtime_error("FitRotoTranslation: zero-norm eigenvector.");
    q_out[0]/=n; q_out[1]/=n; q_out[2]/=n; q_out[3]/=n;
}

inline FitRTResult FitRotoTranslation(const std::vector<Point>& A,
                                      const std::vector<Point>& B)
{
    if(A.size() != B.size())
        throw std::runtime_error("FitRotoTranslation: A and B sizes differ.");
    if(A.empty())
        throw std::runtime_error("FitRotoTranslation: empty point clouds.");

    const std::size_t n = A.size();
    for(std::size_t i=0;i<n;++i){
        if(A[i].coords.size()<3 || B[i].coords.size()<3)
            throw std::runtime_error("FitRotoTranslation: point with <3 coordinates.");
    }

    double ca[3]={0,0,0}, cb[3]={0,0,0};
    for(std::size_t i=0;i<n;++i){
        ca[0]+=A[i].coords[0]; ca[1]+=A[i].coords[1]; ca[2]+=A[i].coords[2];
        cb[0]+=B[i].coords[0]; cb[1]+=B[i].coords[1]; cb[2]+=B[i].coords[2];
    }
    ca[0]/=n; ca[1]/=n; ca[2]/=n;
    cb[0]/=n; cb[1]/=n; cb[2]/=n;

    double S[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(std::size_t i=0;i<n;++i){
        const double ax=A[i].coords[0]-ca[0];
        const double ay=A[i].coords[1]-ca[1];
        const double az=A[i].coords[2]-ca[2];

        const double bx=B[i].coords[0]-cb[0];
        const double by=B[i].coords[1]-cb[1];
        const double bz=B[i].coords[2]-cb[2];

        S[0][0]+=ax*bx; S[0][1]+=ax*by; S[0][2]+=ax*bz;
        S[1][0]+=ay*bx; S[1][1]+=ay*by; S[1][2]+=ay*bz;
        S[2][0]+=az*bx; S[2][1]+=az*by; S[2][2]+=az*bz;
    }

    const double Sxx=S[0][0], Sxy=S[0][1], Sxz=S[0][2];
    const double Syx=S[1][0], Syy=S[1][1], Syz=S[1][2];
    const double Szx=S[2][0], Szy=S[2][1], Szz=S[2][2];
    const double trace = Sxx + Syy + Szz;

    double N[4][4];
    N[0][0]= trace;      N[0][1]= Syz - Szy;  N[0][2]= Szx - Sxz;  N[0][3]= Sxy - Syx;
    N[1][0]= N[0][1];    N[1][1]= Sxx - Syy - Szz;  N[1][2]= Sxy + Syx;  N[1][3]= Szx + Sxz;
    N[2][0]= N[0][2];    N[2][1]= N[1][2];          N[2][2]= -Sxx + Syy - Szz; N[2][3]= Syz + Szy;
    N[3][0]= N[0][3];    N[3][1]= N[1][3];          N[3][2]= N[2][3];          N[3][3]= -Sxx - Syy + Szz;

    double q[4];
    LargestEigenvector4x4_Symmetric_Jacobi(N, q);

    FitRTResult out;
    out.n = n;

    QuaternionToR(q, out.T.R);

    // t = cb - R*ca
    const double Rca0 = out.T.R[0][0]*ca[0] + out.T.R[0][1]*ca[1] + out.T.R[0][2]*ca[2];
    const double Rca1 = out.T.R[1][0]*ca[0] + out.T.R[1][1]*ca[1] + out.T.R[1][2]*ca[2];
    const double Rca2 = out.T.R[2][0]*ca[0] + out.T.R[2][1]*ca[1] + out.T.R[2][2]*ca[2];

    out.T.t[0] = cb[0] - Rca0;
    out.T.t[1] = cb[1] - Rca1;
    out.T.t[2] = cb[2] - Rca2;

    // RMS residual using authoritative R,t
    double ss=0.0;
    for(std::size_t i=0;i<n;++i){
        double xo,yo,zo;
        out.T.apply(A[i].coords[0], A[i].coords[1], A[i].coords[2], xo,yo,zo);
        const double dx = xo - B[i].coords[0];
        const double dy = yo - B[i].coords[1];
        const double dz = zo - B[i].coords[2];
        ss += dx*dx + dy*dy + dz*dz;
    }
    out.rms = std::sqrt(ss / (3.0*n));

    return out;
}

#endif // FIT_ROTOTRANSLATION_H

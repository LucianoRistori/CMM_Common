#include "MapWarm2Cold.h"
#include <cmath>

std::pair<double,double>
MapWarm2Cold(double xw, double yw,
             double tx, double ty,
             double s, double theta)
{
    const double c = std::cos(theta);
    const double si = std::sin(theta);

    const double xr =  c * xw - si * yw;
    const double yr =  si * xw + c * yw;

    return {
        tx + s * xr,
        ty + s * yr
    };
}

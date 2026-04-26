#include "InitialValues.h"
#include <random>
#include <cmath>

template <int dim, int n>
double InitialValues<dim, n>::value(const Point<dim> &p,
                                    const unsigned int component) const
{
 static std::mt19937 gen(std::random_device{}());
static std::normal_distribution<double> dist(0.0, 1.0);

  double noise = dist(gen) * 1e-6;

const double r = p.norm();
const double theta = std::atan2(p[1], p[0]);

if (r < 5.0)
{
    if (component == 0)
    {
        if (theta >= -M_PI/3.0 && theta < M_PI/3.0)
            return 0.8 + noise;
        else
            return 0.0 + noise;
    }

    else if (component == 1)
    {
        if (theta >= M_PI/3.0 && theta < M_PI)
            return 0.8 + noise;
        else
            return 0.0 + noise;
    }

    else if (component == 2)
    {
        if (theta >= -M_PI && theta < -M_PI/3.0)
            return 0.8 + noise;
        else
            return 0.0 + noise;
    }
}

return 0.0;
}

// !!! DŮLEŽITÉ: Explicitní instanciace !!!
template class InitialValues<2, 3>;
// Přidej další podle potřeby v model.cc 
 
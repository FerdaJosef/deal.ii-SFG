#include "InitialValues.h"
#include <random>
#include <cmath>

template <int dim, int n>
double InitialValues<dim, n>::value(const Point<dim> &p,
                                    const unsigned int component) const
{
  (void)p;
  (void)component;
  static std::mt19937 gen(std::random_device{}());
  static std::normal_distribution<double> dist(0.0, 1.0);

  double noise = dist(gen) * 1e-6;

  const double r = p.norm();
  const double theta = std::atan2(p[1], p[0]);

  if (r < 7.5)
  {
      if (component == 0)
      {
        return 0.99 + noise;
      }

      else if (component == 1)
      {
        return 0.0;
      }

      if (component == 2)
      {
        return 0.0;
      }

      else if (component == 3)
      {
        return 0.0;
      }
  }
    return 0.0;
}

// !!! DŮLEŽITÉ: Explicitní instanciace !!!
template class InitialValues<2, 4>;
// Přidej další podle potřeby v model.cc 
 
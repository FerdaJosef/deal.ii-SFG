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

  return 0.0*noise;
}

// !!! DŮLEŽITÉ: Explicitní instanciace !!!
template class InitialValues<3, 2>;
// Přidej další podle potřeby v model.cc 
 
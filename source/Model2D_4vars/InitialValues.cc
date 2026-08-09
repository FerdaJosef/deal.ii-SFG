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

  // Increased noise scale to match the rough interface on the left
  const double noise = dist(gen) * 0.05;
  const double r = p.norm();

  if (r < 5.0)
  {
    const double theta = std::atan2(p[1], p[0]);

    // Component 0: Red -> Top-left sector [90°, 210°]
    if (component == 0)
    {
      if (theta >= M_PI / 2.0 || theta < -5.0 * M_PI / 6.0)
        return 0.99 + noise;
      else
        return 0.0;
    }
    // Component 1: Green -> Bottom sector [-150°, -30°]
    else if (component == 1)
    {
      if (theta >= -5.0 * M_PI / 6.0 && theta < -M_PI / 6.0)
        return 0.99 + noise;
      else
        return 0.0;
    }
    // Component 2: Blue -> Top-right sector [-30°, 90°]
    else if (component == 2)
    {
      if (theta >= -M_PI / 6.0 && theta < M_PI / 2.0)
        return 0.99 + noise;
      else
        return 0.0;
    }
  }

  // Matrix background retains small noise field
  return noise;
}

template class InitialValues<2, 4>;
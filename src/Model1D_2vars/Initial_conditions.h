#pragma once

#include <deal.II/base/tensor.h>
#include <vector>
#include <random>

template <int dim, int n>
class InitialCondition : public Function<dim>
{
  public:
  InitialCondition(const double time = 0.)
  : Function<dim>(n, time)
  {}
  
virtual double value(const Point<dim> &p,
                     const unsigned int component = 0) const override
{
  static std::mt19937 gen(std::random_device{}());
  static std::normal_distribution<double> dist(0.0, 1.0);

  const double noise = dist(gen) * 2e-2;

  if (p[0] <= 10.0)
  {
    if (component == 0) {
      double hodnota = 0.5 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }


    else if (component == 1)
    {
      double hodnota =  0.5 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

  }

  else if (p[0] <= 17.0)
  {
    double hodnota =  0.0 + noise;
    return std::max(0.0, std::min(1.0, hodnota));
  }

  else if (p[0] <= 25.0)
  {
    if (component == 0) {
      double hodnota = 0.0 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }


    else if (component == 1) {
      double hodnota =  0.5 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

  }

  else if (p[0] <= 33.0)
  {
    if (component == 0) {
      double hodnota =  0.5 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }


    else if (component == 1) {
      double hodnota = 0.0 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

  }

  else if (p[0] <= 37.0)
  {
    double hodnota = 0.0 + noise;
    return std::max(0.0, std::min(1.0, hodnota));
  }

  else if (p[0] <= 40.0)
  {
    if (component == 0) {
      double hodnota =  1.0 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

    else if (component == 1) {
      double hodnota =  0.0 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

  }

  else if (p[0] <= 43.0)
  {
    if (component == 0) {
      double hodnota =  0.0 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

    else if (component == 1) {
      double hodnota =  1.0 + noise;
      return std::max(0.0, std::min(1.0, hodnota));
    }

  }

  return 0.0;
}
};
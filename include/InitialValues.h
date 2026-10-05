#ifndef INITIAL_VALUES_H
#define INITIAL_VALUES_H

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

template <int dim, int n>
class InitialValues : public Function<dim>
{
public:
  InitialValues(const double time = 0.)
    : Function<dim>(n, time)
  {}

  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
};

#endif
#pragma once

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/tensor.h>
#include <vector>
#include <random>

template <int dim, int n>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide(double delta_t)
    : Function<dim>(n)
    , delta_t(delta_t)
  {}

  virtual double value(const Point<dim>  &p,
                      const unsigned int component) const override;

private:
  double delta_t;
};

template <int dim, int n>
double RightHandSide<dim, n>::value(const Point<dim> &p, const unsigned int component) const
{
    (void)p;
    (void)component;

    static std::mt19937 gen(std::random_device{}());
    static std::normal_distribution<double> dist(0.0,1.0);

    return 0.0;

}

template <int dim, int n>
class RandomRHS
{
public:
  static void generate_rhs(
    double delta_t,
    unsigned int n_q_points,
    const DoFHandler<dim> &dof_handler,
    std::vector<std::vector<Tensor<1,n>>> &rhs_values);
};

template <int dim, int n>
void RandomRHS<dim, n>::generate_rhs(
  double delta_t,
  const unsigned int n_q_points,
  const DoFHandler<dim> &dof_handler,
  std::vector<std::vector<Tensor<1,n>>> &rhs_values)
{
  static thread_local std::mt19937 gen(std::random_device{}());
  static thread_local std::normal_distribution<double> dist(0.0,1.0);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    const unsigned int cell_id = cell->active_cell_index();

    for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int r=0; r < n; r++)
      {
        rhs_values[cell_id][q][r] =
            1e-6*dist(gen)/std::sqrt(delta_t);
      }
    }
  }
}
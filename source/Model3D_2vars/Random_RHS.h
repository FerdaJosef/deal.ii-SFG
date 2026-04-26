#include "model.h"

template <int dim, int n>
void Step3<dim, n>::generate_rhs()
{
  static std::mt19937 gen(std::random_device{}());
  static std::normal_distribution<double> dist(0.0,1.0);

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
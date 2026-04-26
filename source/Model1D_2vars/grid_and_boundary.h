#include "model.h"

template <int dim, int n>
void Step3<dim, n>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, left_lim, right_lim, true);
  triangulation.refine_global(
    n_refinements
  );

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}

template <int dim, int n>
void Step3<dim, n>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  constraints.clear();

  if constexpr (dim == 2)
  {
    DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
    DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
  }

  if constexpr (dim == 3)
  {
    DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
    DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
    DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, constraints);
  }

  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  oldsolution.reinit(dof_handler.n_dofs());
  newton_iterate.reinit(dof_handler.n_dofs());

  rhs_values.resize(triangulation.n_active_cells(),
                  std::vector<Tensor<1,n>>(n_q_points));
}
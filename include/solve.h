#include "model.h"

template <int dim, int n>
double Step3<dim, n>::determine_step_length() const
{
  return 1.0;
}

template <int dim, int n>
void Step3<dim, n>::solve()
{
  SolverControl            solver_control(20000, 1e-6 * system_rhs.l2_norm());
  SolverGMRES<Vector<double>> solver(solver_control);

  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);

  solver.solve(system_matrix, newton_iterate, system_rhs, preconditioner);

  constraints.distribute(newton_iterate);

  const double alpha = determine_step_length();

  solution.add(alpha, newton_iterate);

  constraints.distribute(solution);

  solver_iteration = solver_control.last_step();

  std::cout << solver_iteration
            << " iterations needed to obtain convergence." << std::endl;
}
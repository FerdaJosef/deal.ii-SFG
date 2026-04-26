#include "model.h"

template <int dim, int n>
Step3<dim, n>::Step3(ParameterHandler &param)
  : prm(param)
  , fe(FE_Q<dim>(1), n)
  , dof_handler(triangulation)
  , n_q_points(QGauss<dim>(fe.degree + 1).size())
  , newton_iteration(0)
{
  prm.enter_subsection("Mesh & geometry parameters");

  n_refinements = prm.get_integer("Number of refinements");
  left_lim = prm.get_double("Mesh left limit");
  right_lim = prm.get_double("Mesh right limit");

  prm.leave_subsection();

  prm.enter_subsection("Temporal parameters");
  {
    time           = prm.get_double("Initial time");
    delta_t          = prm.get_double("Initial time step");
    final_time           = prm.get_double("Final time");
    timestep_number           = prm.get_integer("Time step number");
    dt_max         = prm.get_double("Max time step");
    dt_min         = prm.get_double("Max time step");
    max_multiplier = prm.get_double("Max multiplier");
    min_multiplier = prm.get_double("Min multiplier");
    optimal_it = prm.get_integer("Optimal iterations");
    max_it = prm.get_integer("Max iterations");
  }
  prm.leave_subsection();

}

template class Step3<1,2>;
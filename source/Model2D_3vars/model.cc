#include "model.h"

template <int dim, int n>
Step3<dim, n>::Step3()
  : fe(FE_Q<dim>(1), n)
  , dof_handler(triangulation)
  , n_q_points(QGauss<dim>(fe.degree + 1).size())
  , time(0.0)
  , final_time(800.0)
  , delta_t(1e-4)
  , timestep_number(0)
  , max_it(10)
  , max_multiplier(2.0)
  , min_multiplier(0.5)
  , optimal_it(6)
  , dt_max(5.0)
  , dt_min(1e-7)
  , newton_iteration(0)
{}
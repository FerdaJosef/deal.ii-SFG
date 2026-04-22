#include "source/Model2D_2vars/model.h"

template <int dim, int n>

void Step3<dim, n>::run()
{

  make_grid();
  setup_system();

  std::cout << n_q_points << std::endl;

  VectorTools::project(dof_handler,
                      nonzero_constraints,
                      QGauss<dim>(fe.degree + 1),
                      ExactSolution<dim, n>(0.0),
                      solution);

  nonzero_constraints.distribute(solution);

  //std::cout << "The difference is " << L2_error << std::endl;

  output_results();

  std::cout << time << std::endl;
  while (time < final_time)
  {
    make_timestep();

    bool value = time_step_update();
    if (!value)
      continue;
  
    output_results();
  }

  std::cout << solution.linfty_norm() << std::endl;
}
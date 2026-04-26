#include "model.h"

template <int dim, int n>
bool Step3<dim, n>::time_step_update()
{
    if (newton_iteration == max_it)
    {
        std::cout << "Newton failed. Reducing time step." << std::endl;
        solution = oldsolution;

        time -= delta_t;
        timestep_number--;

        delta_t *= 0.5;
        if (delta_t < dt_min)
            throw std::runtime_error("Time step too small!");

        return false;
    }

    if (solver_iteration > 250) {
        std::cout << "GMRES failed. Reducing time step." << std::endl;
        solution = oldsolution;

        time -= delta_t;
        timestep_number--;

        delta_t *= 0.5;
        if (delta_t < dt_min)
            throw std::runtime_error("Solver failed!");

        return false;
    }

    if (newton_iteration <= optimal_it) {
      // self.dt.assign(min(dt_max, float(self.dt) * min(1.6, optimal_it / (iterations + 0.001))));
        delta_t = std::min(dt_max, delta_t*std::min(max_multiplier, double(optimal_it/(newton_iteration + 0.001))));
    }
    else {
      // self.dt.assign(float(self.dt) * max(0.5, optimal_it / iterations));
      delta_t = delta_t*std::max(min_multiplier, double(optimal_it / newton_iteration));
    }

    if (time + delta_t > final_time) {
        delta_t = final_time - time;
    }
      
    else if (delta_t < 1e-6) {
        throw std::invalid_argument("Time step too small!");
    }

    std::cout << newton_iteration << std::endl;
    
    std::cout << "Our time step is " << delta_t << std::endl;

    return true;
}

template<int dim, int n>
void Step3<dim, n>::make_timestep()
{
    time+=delta_t;
    timestep_number++;

    generate_rhs();

    std::cout<< "Time: " << time << std::endl;
    oldsolution = solution;
    newton_iteration = 0;
    double residual_norm = 1.0;
    while (residual_norm > 1e-10 && newton_iteration < max_it)
    {
      //newton_iterate = 0.0;
      assemble_system();
      residual_norm = system_rhs.l2_norm();

      std::cout << "The norm of our solution is: " << residual_norm << std::endl;

      solve();

      //std::cout << "  Residual: " << residual_norm << std::endl;
      newton_iteration++;
    }
}
#include "model.h"
#include "InitialValues.h"
#include "RandomField.h"

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

};

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

  random_field.reinit(triangulation.n_active_cells(), n_q_points);
};

template <int dim, int n>
Step3<dim, n>::AssemblyScratchData::AssemblyScratchData(const FiniteElement<dim> &fe)
  :
  fe_values(fe,
            QGauss<dim>(fe.degree + 1),
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values)
{
  const unsigned int n_q_points =
      fe_values.get_quadrature().size();

  values_newton.resize(n_q_points, Vector<double>(n));
  values_old.resize(n_q_points, Vector<double>(n));

  gradients_newton.resize(n_q_points,
                          std::vector<Tensor<1,dim>>(n));

  dPsiDu.reinit(n);
  dPsiDu2.reinit(n, n);

  dPsidGradU.resize(n);

  dPsidUdGradU.resize(n, std::vector<Tensor<1,dim>>(n));

  dPsidGradU2.resize(n, std::vector<Tensor<2,dim>>(n));

  acegen_scratch.resize(256);
}

template <int dim, int n>
Step3<dim, n>::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
    , values_newton(scratch_data.values_newton)
    , values_old(scratch_data.values_old)
    , gradients_newton(scratch_data.gradients_newton)
    , acegen_scratch(scratch_data.acegen_scratch)
    , dPsiDu(scratch_data.dPsiDu)
    , dPsidGradU(scratch_data.dPsidGradU)
    , dPsiDu2(scratch_data.dPsiDu2)
    , dPsidUdGradU(scratch_data.dPsidUdGradU)
    , dPsidGradU2(scratch_data.dPsidGradU2)
  {}

  template <int dim, int n>
void Step3<dim, n>::assemble_system()
{
  system_matrix = 0;
  system_rhs = 0;

  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &Step3<dim, n>::local_assemble_system,
                  &Step3<dim, n>::copy_local_to_global,
                  AssemblyScratchData(fe),
                  AssemblyCopyData());
}

template <int dim, int n>
void Step3<dim, n>::local_assemble_system(
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        AssemblyScratchData                                  &scratch_data,
        AssemblyCopyData                                     &copy_data)
{
  const QGauss<dim> quadrature_formula(fe.degree + 1);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  //n_q_points    = quadrature_formula.size();

  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);

  copy_data.local_dof_indices.resize(dofs_per_cell);
 
  scratch_data.fe_values.reinit(cell);

  //======= ACEGEN=======

  scratch_data.fe_values.reinit(cell);

  copy_data.cell_matrix = 0;
  copy_data.cell_rhs    = 0;

  scratch_data.fe_values.get_function_gradients(solution, scratch_data.gradients_newton);
  scratch_data.fe_values.get_function_values(oldsolution, scratch_data.values_old);
  scratch_data.fe_values.get_function_values(solution, scratch_data.values_newton);

  const unsigned int cell_id = cell->active_cell_index();

  //right_hand_side(fe_values.get_quadrature_points(), rhs_values);

  for (const unsigned int q_index : scratch_data.fe_values.quadrature_point_indices())
    {
    //const auto &x_q = fe_values.quadrature_point(q_index);

    equation<dim, n>(scratch_data.acegen_scratch,
            scratch_data.values_newton[q_index],
            scratch_data.values_old[q_index],
        scratch_data.gradients_newton[q_index],
        scratch_data.dPsiDu,
        scratch_data.dPsidGradU,
        scratch_data.dPsiDu2,
        scratch_data.dPsidUdGradU,
        scratch_data.dPsidGradU2, // 2, 2, 3, 3
        &this->delta_t 
      );

      const auto &rhs_val = random_field.get_value(cell->active_cell_index(), q_index);

      const auto &sd = scratch_data;
      for (const unsigned int i : scratch_data.fe_values.dof_indices())
        {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        
        // number of degrees of freedom is equql to number "geometric dof" = support points
        // function fe.system_to_component_index(i) returns a pair:
        // first is the corresponding component of vector system
        // second is the index of support point
        // shape value known about it through index i
        // ted tim vzdycky zredukuju dimenzi objektu a prevedu to na skalarni pripad

          copy_data.cell_rhs(i) -= (sd.fe_values.shape_value(i, q_index) * 
                          sd.dPsiDu[component_i]

                      + sd.fe_values.shape_grad(i, q_index) * sd.dPsidGradU[component_i]
                          + sd.fe_values.shape_value(i,q_index) * rhs_val[component_i]) *
                          sd.fe_values.JxW(q_index);
          }
      for (const unsigned int i : scratch_data.fe_values.dof_indices())
      {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        
          for (const unsigned int j : scratch_data.fe_values.dof_indices())
          {
            const unsigned int component_j = fe.system_to_component_index(j).first;

            copy_data.cell_matrix(i, j) +=
                (sd.fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    sd.dPsidGradU2[component_i][component_j] *           // dPsi/d2(grad u)
                sd.fe_values.shape_grad(j, q_index)    // grad phi_j(x_q)
                + 
                sd.fe_values.shape_value(i, q_index) * // phi_i(x_q)
                    sd.dPsiDu2(component_i, component_j) *           // dPsi/d2( u)
                sd.fe_values.shape_value(j, q_index)     // phi_j(x_q)
                +
                sd.fe_values.shape_value(i, q_index) * //  phi_i(x_q)
                    sd.dPsidUdGradU[component_i][component_j] *           // dPsi/d(grad u)du
                sd.fe_values.shape_grad(j, q_index)     // grad phi_j(x_q)
                +
                sd.fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    sd.dPsidUdGradU[component_i][component_j] *           // dPsi/d(grad u)du
                sd.fe_values.shape_value(j, q_index)     // phi_j(x_q)
                ) *
                sd.fe_values.JxW(q_index);           // dx
        }
      }
    }
  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim, int n>
void Step3<dim, n>::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(
      copy_data.cell_matrix,
      copy_data.cell_rhs,
      copy_data.local_dof_indices,
      system_matrix,
      system_rhs);
  }

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
    random_field.generate(dof_handler.n_dofs(), delta_t, 1e-5, triangulation.n_active_cells());

    time+=delta_t;
    timestep_number++;

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

template <int dim, int n>
void Step3<dim, n>::output_results() const
{
  static std::vector<std::pair<double, std::string>> times_and_names;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  prm.enter_subsection("Output parameters");
 
  const std::string output_filename = prm.get("Output filename");
  data_out.parse_parameters(prm);
 
  prm.leave_subsection();

  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  const std::string filename = output_filename + Utilities::int_to_string(timestep_number) + ".vtu";
  std::ofstream output(filename);
  data_out.write(output, DataOutBase::vtu);

  std::cout << "Output written to " << filename << std::endl;
  times_and_names.push_back(
                {time, filename});

  std::ofstream pvd_output (output_filename);
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}

template <int dim, int n>

void Step3<dim, n>::run()
{

  make_grid();
  setup_system();

  std::cout << n_q_points << std::endl;

  VectorTools::project(dof_handler,
                      constraints,
                      QGauss<dim>(fe.degree + 1),
                      InitialValues<dim, n>(0.0),
                      solution);

  constraints.distribute(solution);

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

template class Step3<3,2>;
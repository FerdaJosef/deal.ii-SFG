/**
 * @file model.cc
 * @brief Implementation of Step3 for Stochastic Phase-Field modeling using MPI.
 */

#include "model.h"
#include "InitialValues.h"

template <int dim, int n>
Step3<dim, n>::Step3(ParameterHandler &param)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, (this_mpi_process == 0))
  , prm(param)
  , computing_timer(mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times)
  , fe(FE_Q<dim>(1), n)
  , dof_handler(triangulation)
  , n_q_points(QGauss<dim>(fe.degree + 1).size())
  , newton_iteration(0)
{
  prm.enter_subsection("Mesh & geometry parameters");
  n_refinements = prm.get_integer("Number of refinements");
  left_lim      = prm.get_double("Mesh left limit");
  right_lim     = prm.get_double("Mesh right limit");
  prm.leave_subsection();

  prm.enter_subsection("Temporal parameters");
  {
    time            = prm.get_double("Initial time");
    delta_t         = prm.get_double("Initial time step");
    final_time      = prm.get_double("Final time");
    timestep_number = prm.get_integer("Time step number");
    dt_max          = prm.get_double("Max time step");
    dt_min          = prm.get_double("Min time step");
    max_multiplier  = prm.get_double("Max multiplier");
    min_multiplier  = prm.get_double("Min multiplier");
    optimal_it      = prm.get_integer("Optimal iterations");
    max_it          = prm.get_integer("Max iterations");
  }
  prm.leave_subsection();

  prm.enter_subsection("Solver");
  {
    linear_residual      = prm.get_double("Linear system error");
    max_linear_iteration = prm.get_integer("Max linear solver iterations");
  }
  prm.leave_subsection();
}

template <int dim, int n>
void Step3<dim, n>::make_grid()
{
  TimerOutput::Scope scope(computing_timer, "Making grid");
  GridGenerator::hyper_cube(triangulation, left_lim, right_lim, true);
  triangulation.refine_global(n_refinements);

  pcout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}

template <int dim, int n>
void Step3<dim, n>::setup_system()
{
  TimerOutput::Scope timing_section(computing_timer, "Setting up our system");

#ifdef DEAL_II_WITH_METIS
  GridTools::partition_triangulation(n_mpi_processes, triangulation);
#else
  GridTools::partition_triangulation_zorder(n_mpi_processes, triangulation);
#endif

  dof_handler.distribute_dofs(fe);
  DoFRenumbering::subdomain_wise(dof_handler);

  const std::vector<IndexSet> locally_owned_dofs_per_proc =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
  locally_owned_dofs = locally_owned_dofs_per_proc[this_mpi_process];
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

  constraints.clear();
  // FIXED: Updated to use the non-deprecated two-argument reinit
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

  if constexpr (dim == 2)
  {
    DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
    DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
  }
  else if constexpr (dim == 3)
  {
    DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
    DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
    DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, constraints);
  }
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp, locally_owned_dofs, mpi_communicator, locally_relevant_dofs);

  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

  // Non-ghosted vectors (for calculations and linear algebra)
  distributed_solution.reinit(locally_owned_dofs, mpi_communicator);
  distributed_old_solution.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  newton_iterate.reinit(locally_owned_dofs, mpi_communicator);

  // Ghosted vectors (read-only for cell evaluation & output)
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  oldsolution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  random_field.reinit(triangulation.n_active_cells(), n_q_points);
}

template <int dim, int n>
void Step3<dim, n>::assemble_system()
{
  TimerOutput::Scope timing_section(computing_timer, "Assembling");

  system_matrix = 0;
  system_rhs    = 0;

  FEValues<dim> fe_values(fe,
                          QGauss<dim>(fe.degree + 1),
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Quadrature point evaluation buffers
  std::vector<Vector<double>> local_values_newton(n_q_points, Vector<double>(n));
  std::vector<Vector<double>> local_values_old(n_q_points, Vector<double>(n));
  std::vector<std::vector<Tensor<1, dim>>> local_gradients_newton(
    n_q_points, std::vector<Tensor<1, dim>>(n));

  // AceGen output containers
  Vector<double>                          dPsiDu(n);
  FullMatrix<double>                      dPsiDu2(n, n);
  std::vector<Tensor<1, dim>>             dPsidGradU(n);
  std::vector<std::vector<Tensor<1, dim>>> dPsidUdGradU(n, std::vector<Tensor<1, dim>>(n));
  std::vector<std::vector<Tensor<2, dim>>> dPsidGradU2(n, std::vector<Tensor<2, dim>>(n));
  std::vector<double>                     acegen_scratch(256);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->subdomain_id() == this_mpi_process)
    {
      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.reinit(cell);
      fe_values.get_function_values(solution, local_values_newton);
      fe_values.get_function_values(oldsolution, local_values_old);
      fe_values.get_function_gradients(solution, local_gradients_newton);

      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
      {
        equation<dim, n>(acegen_scratch,
                         local_values_newton[q_index],
                         local_values_old[q_index],
                         local_gradients_newton[q_index],
                         dPsiDu,
                         dPsidGradU,
                         dPsiDu2,
                         dPsidUdGradU,
                         dPsidGradU2,
                         &this->delta_t);

        const auto &rhs_val = random_field.get_value(cell->active_cell_index(), q_index);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe.system_to_component_index(i).first;

          cell_rhs(i) -= (fe_values.shape_value(i, q_index) * dPsiDu[component_i] +
                          fe_values.shape_grad(i, q_index) * dPsidGradU[component_i] +
                          fe_values.shape_value(i, q_index) * rhs_val[component_i]) *
                         fe_values.JxW(q_index);

          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            const unsigned int component_j = fe.system_to_component_index(j).first;

            cell_matrix(i, j) +=
              (fe_values.shape_grad(i, q_index) * dPsidGradU2[component_i][component_j] *
                 fe_values.shape_grad(j, q_index) +
               fe_values.shape_value(i, q_index) * dPsiDu2(component_i, component_j) *
                 fe_values.shape_value(j, q_index) +
               fe_values.shape_value(i, q_index) * dPsidUdGradU[component_i][component_j] *
                 fe_values.shape_grad(j, q_index) +
               fe_values.shape_grad(i, q_index) * dPsidUdGradU[component_j][component_i] *
                 fe_values.shape_value(j, q_index)) *
              fe_values.JxW(q_index);
          }
        }
      }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             cell_rhs,
                                             local_dof_indices,
                                             system_matrix,
                                             system_rhs);
    }
  }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <int dim, int n>
bool Step3<dim, n>::time_step_update()
{
  if (newton_iteration == max_it)
  {
    pcout << "Newton failed. Reducing time step." << std::endl;
    distributed_solution = distributed_old_solution;
    solution = distributed_solution;

    time -= delta_t;
    timestep_number--;

    delta_t *= 0.5;
    if (delta_t < dt_min)
      throw std::runtime_error("Time step too small!");

    return false;
  }

  if (solver_iteration > max_linear_iteration)
  {
    pcout << "Linear solver failed. Reducing time step." << std::endl;
    distributed_solution = distributed_old_solution;
    solution = distributed_solution;

    time -= delta_t;
    timestep_number--;

    delta_t *= 0.5;

    if (delta_t < dt_min)
      throw std::runtime_error("Solver failed!");

    return false;
  }

  if (newton_iteration <= optimal_it)
  {
    delta_t = std::min(dt_max, delta_t * std::min(max_multiplier, double(optimal_it / (newton_iteration + 0.001))));
  }
  else
  {
    delta_t = delta_t * std::max(min_multiplier, double(optimal_it / newton_iteration));
  }

  if (time + delta_t > final_time)
  {
    delta_t = final_time - time;
  }
  else if (delta_t < 1e-6)
  {
    throw std::invalid_argument("Time step too small!");
  }

  pcout << "Newton iterations: " << newton_iteration << std::endl;
  pcout << "Our time step is " << delta_t << std::endl;

  return true;
}

template <int dim, int n>
void Step3<dim, n>::sync_solution_and_assemble(const PETScWrappers::MPI::Vector &u)
{
  // u is KINSOL's (non-ghosted) evaluation point; push it into our
  // distributed vector, then refresh the ghosted copy assemble_system() reads from.
  distributed_solution = u;
  solution = distributed_solution;
  assemble_system();
}

template <int dim, int n>
void Step3<dim, n>::make_timestep()
{
  random_field.generate(triangulation.n_active_cells(), n_q_points, delta_t, 1e-6);

  time += delta_t;
  timestep_number++;

  pcout << "Time: " << time << std::endl;

  distributed_old_solution = distributed_solution;
  oldsolution = distributed_old_solution;

  newton_iteration = 0;

  SUNDIALS::KINSOL<PETScWrappers::MPI::Vector>::AdditionalData additional_data;
  // 1. Relax tolerances to realistic values for 264k DOFs
  additional_data.function_tolerance            = 1e-8; 
  additional_data.step_tolerance                = 1e-8;
  additional_data.maximum_non_linear_iterations = max_it;
  
  // 2. Disable line-search backtracking to prevent endless residual loops
  additional_data.strategy = SUNDIALS::KINSOL<PETScWrappers::MPI::Vector>::AdditionalData::newton;

  SUNDIALS::KINSOL<PETScWrappers::MPI::Vector> nonlinear_solver(additional_data, mpi_communicator);

  nonlinear_solver.reinit_vector = [&](PETScWrappers::MPI::Vector &v)
  {
    v.reinit(locally_owned_dofs, mpi_communicator);
  };

nonlinear_solver.residual =
    [&](const PETScWrappers::MPI::Vector &u, PETScWrappers::MPI::Vector &F) -> int
  {
    sync_solution_and_assemble(u);
    F = system_rhs;
    F *= -1.0; 
    ++newton_iteration;
    
    double res_norm = F.l2_norm();
    pcout << "    Newton Iteration " << newton_iteration 
          << " | Residual ||F||_2: " << res_norm << std::endl;
          
    if (std::isnan(res_norm)) {
        pcout << "\n[!] NaN detected in Residual Assembly! The material model evaluated out-of-bounds." << std::endl;
    }
          
    return 0;
  };

  nonlinear_solver.solve_with_jacobian =
    [&](const PETScWrappers::MPI::Vector &rhs, PETScWrappers::MPI::Vector &dst,
        const double /*tolerance*/) -> int
  {
    SolverControl solver_control(max_linear_iteration, std::max(linear_residual * rhs.l2_norm(), 1e-10));

    PETScWrappers::SparseDirectMUMPS direct_solver(solver_control);

    pcout << "Frobenius: " << system_matrix.frobenius_norm() << std::endl;
    pcout << "RHS norm: " << rhs.l2_norm() << std::endl;

    direct_solver.solve(system_matrix, dst, rhs);
    
    constraints.distribute(dst);

    double step_norm = dst.l2_norm();
    pcout << "      -> Linear Solve step ||Delta u||_2: " << step_norm << std::endl;
    
    if (std::isnan(step_norm)) {
        pcout << "\n[!] NaN detected in Linear Solve! The Jacobian matrix is likely singular." << std::endl;
    }

    solver_iteration = solver_control.last_step();
    return 0;
  };

  nonlinear_solver.setup_jacobian =
    [&](const PETScWrappers::MPI::Vector &u, const PETScWrappers::MPI::Vector & /*F*/) -> int
  {
    (void)u;
    return 0;
  };

  try
  {
    nonlinear_solver.solve(distributed_solution);
    solution = distributed_solution;
    pcout << "  -> KINSOL converged successfully in " << newton_iteration << " iterations." << std::endl;
  }
  catch (const std::exception &e)
  {
    pcout << "\n[!] KINSOL failed to converge: " << e.what() << std::endl;
    pcout << "[!] Redirecting to adaptive time-step reduction logic..." << std::endl;
    newton_iteration = max_it;
  }
}

template <int dim, int n>
void Step3<dim, n>::output_results() const
{
  TimerOutput::Scope timing_section(computing_timer, "Outputting");

  const Vector<double> localized_solution(solution);

  if (this_mpi_process == 0)
  {
    static std::vector<std::pair<double, std::string>> times_and_names;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    prm.enter_subsection("Output parameters");
    const std::string output_filename = prm.get("Output filename");
    data_out.parse_parameters(prm);
    prm.leave_subsection();

    data_out.add_data_vector(localized_solution, "solution");

    std::vector<unsigned int> partition_int(triangulation.n_active_cells());
    GridTools::get_subdomain_association(triangulation, partition_int);

    const Vector<double> partitioning(partition_int.begin(),
                                      partition_int.end());

    data_out.add_data_vector(partitioning, "partitioning");

    data_out.build_patches();

    const std::string filename = output_filename + Utilities::int_to_string(timestep_number) + ".vtu";
    const std::string vtu_basename = std::filesystem::path(filename).filename().string();

    std::ofstream output(filename);
    data_out.write(output, DataOutBase::vtu);

    pcout << "Output written to " << filename << std::endl;
    times_and_names.push_back({time, vtu_basename});

    std::ofstream pvd_output(output_filename + ".pvd");
    DataOutBase::write_pvd_record(pvd_output, times_and_names);
  }
}

template <int dim, int n>
void Step3<dim, n>::run()
{
  make_grid();
  setup_system();

  pcout << "Quadrature points per cell: " << n_q_points << std::endl;

  VectorTools::project(dof_handler,
                       constraints,
                       QGauss<dim>(fe.degree + 1),
                       InitialValues<dim, n>(0.0),
                       distributed_solution);

  constraints.distribute(distributed_solution);
  solution = distributed_solution; // Syncs non-ghosted to ghosted vector

  output_results();

  pcout << "Starting time: " << time << std::endl;
  while (time < final_time)
  {
    make_timestep();

    bool value = time_step_update();
    if (!value)
      continue;

    output_results();
  }

  pcout << "L-infinity norm: " << solution.linfty_norm() << std::endl;
}

template class Step3<2, 4>;
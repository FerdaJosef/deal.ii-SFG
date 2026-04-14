/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 1999 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */



#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include "AceGen/2D_3omega/equation_2_3.h"
#include "wrapper.h"

#include <math.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/grid_tools.h>
 
#include <deal.II/dofs/dof_renumbering.h>
 
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
 
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

#include <random>

using namespace dealii;

template <int dim, int n>
class Step3
{
public:
  Step3();

  void run();


private:
  void make_grid();
  void setup_system();

struct AssemblyScratchData
{
  AssemblyScratchData(const FiniteElement<dim> &fe);
  AssemblyScratchData(const AssemblyScratchData &scratch_data);

  FEValues<dim> fe_values;

  std::vector<Vector<double>> values_newton;
  std::vector<Vector<double>> values_old;

  std::vector<std::vector<Tensor<1,dim>>> gradients_newton;

  std::vector<double> acegen_scratch;

  // ===== deal.II (readable layer) =====
  Vector<double> dPsiDu;

  std::vector<Tensor<1,dim>> dPsidGradU;

  FullMatrix<double> dPsiDu2;

  std::vector<std::vector<Tensor<1,dim>>> dPsidUdGradU;

  std::vector<std::vector<Tensor<2,dim>>> dPsidGradU2;
};
 
  struct AssemblyCopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };
  
  void assemble_system();

  void local_assemble_system(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      AssemblyScratchData                                  &scratch,
      AssemblyCopyData                                     &copy_data);
  void copy_local_to_global(const AssemblyCopyData &copy_data);

  void solve();
  void time_step_update();
  double determine_step_length() const;
  void output_results() const;
  void generate_rhs();
  double compute_residual();

  Triangulation<dim> triangulation;
  const FESystem<dim>    fe;
  DoFHandler<dim>    dof_handler;

  AffineConstraints<double> zero_constraints;
  AffineConstraints<double> nonzero_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> oldsolution;
  Vector<double> newton_iterate;
  Vector<double> solution;
  Vector<double> system_rhs;

  const unsigned int n_q_points;
  std::vector<std::vector<Tensor<1,n>>> rhs_values;

  double time;
  double final_time;
  double delta_t;
  unsigned int timestep_number;

  int max_it;
  double max_multiplier;
  double min_multiplier;
  int optimal_it;
  double dt_max; double dt_min;
  int newton_iteration;
  int solver_iteration;
};

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

template <int dim, int n>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide(double time = 0.0)
    : Function<dim>(n, time)
    , period(0.2)
  {}

  virtual double value(const Point<dim>  &p,
                       const unsigned int component) const override;
  
  private:
    const double period;
};


template <int dim, int n>
double RightHandSide<dim, n>::value(const Point<dim> &p, const unsigned int component) const
{
    (void)p;
    (void)component;
    //const double t = this->get_time();

    static std::mt19937 gen(std::random_device{}());
    static std::normal_distribution<double> dist(0.0,1.0);

    //double sample;
    //sample = d(gen);

    if (component == 0)
      return std::sin(M_PI*p[0])*std::sin(M_PI*p[1]);

    else if (component == 1)
      return std::sin(M_PI*p[0])*std::sin(M_PI*p[1]);

  else
    return 0.0;

}

template <int dim, int n>
class ExactSolution : public Function<dim>
{
  public:
  ExactSolution(const double time = 0.)
  : Function<dim>(n, time)
  {}
  
virtual double value(const Point<dim> &p,
                     const unsigned int component = 0) const override
{
  static std::mt19937 gen(std::random_device{}());
  static std::normal_distribution<double> dist(0.0, 1.0);

  double noise = dist(gen) * 1e-6;

const double r = p.norm();
const double theta = std::atan2(p[1], p[0]);

if (r < 5.0)
{
    if (component == 0)
    {
        if (theta >= -M_PI/3.0 && theta < M_PI/3.0)
            return 0.8 + noise;
        else
            return 0.0 + noise;
    }

    else if (component == 1)
    {
        if (theta >= M_PI/3.0 && theta < M_PI)
            return 0.8 + noise;
        else
            return 0.0 + noise;
    }

    else if (component == 2)
    {
        if (theta >= -M_PI && theta < -M_PI/3.0)
            return 0.8 + noise;
        else
            return 0.0 + noise;
    }
}

return 0.0;
}
};

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

template <int dim, int n>
void Step3<dim, n>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -10, 10, true);
  triangulation.refine_global(
    7
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

  zero_constraints.clear();
  VectorTools::interpolate_boundary_values(dof_handler,
                                            0,
                                            Functions::ZeroFunction<dim>(n),
                                            zero_constraints);
  zero_constraints.close();

  nonzero_constraints.clear();

  DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, nonzero_constraints); // x-direction
  DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, nonzero_constraints); // y-direction
  //DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, nonzero_constraints); // z-direction

  nonzero_constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints, true);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  oldsolution.reinit(dof_handler.n_dofs());
  newton_iterate.reinit(dof_handler.n_dofs());

  rhs_values.resize(triangulation.n_active_cells(),
                  std::vector<Tensor<1,n>>(n_q_points));
}

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

  //RightHandSide rhs_function;
  //rhs_function.set_time(this->time);

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
                          + sd.fe_values.shape_value(i,q_index) * rhs_values[cell_id][q_index][component_i]) *
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
    nonzero_constraints.distribute_local_to_global(
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

  nonzero_constraints.distribute(newton_iterate);

  const double alpha = determine_step_length();

  solution.add(alpha, newton_iterate);

  nonzero_constraints.distribute(solution);

  solver_iteration = solver_control.last_step();

  std::cout << solver_iteration
            << " iterations needed to obtain convergence." << std::endl;
}

template <int dim, int n>
void Step3<dim, n>::time_step_update()
{
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
    std::cout << "Our time step is " << delta_t << std::endl;
}

template <int dim, int n>
void Step3<dim, n>::output_results() const
{
  static std::vector<std::pair<double, std::string>> times_and_names;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  const std::string filename = "results/solution-" + Utilities::int_to_string(timestep_number) + ".vtu";
  std::ofstream output(filename);
  data_out.write(output, DataOutBase::vtu);

  std::cout << "Output written to " << filename << std::endl;
  times_and_names.push_back(
                {time, filename});

  std::ofstream pvd_output ("solution1.pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}


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

    if (newton_iteration == max_it)
    {
        std::cout << "Newton failed. Reducing time step." << std::endl;
        solution = oldsolution;

        time -= delta_t;
        timestep_number--;

        delta_t *= 0.5;
        if (delta_t < dt_min)
            throw std::runtime_error("Time step too small!");

        continue;
    }

    if (solver_iteration > 250) {
            std::cout << "GMRES failed. Reducing time step." << std::endl;
        solution = oldsolution;

        time -= delta_t;
        timestep_number--;

        delta_t *= 0.5;
        if (delta_t < dt_min)
            throw std::runtime_error("Solver failed!");

        continue;
    }

    std::cout << newton_iteration << std::endl;
    time_step_update();
    output_results();
  }

  std::cout << solution.linfty_norm() << std::endl;
}



int main()
{
    using namespace dealii;
    try
      {
        MultithreadInfo::set_thread_limit();
  
        Step3<2,3> double_ditch;
        double_ditch.run();
      }
    catch (std::exception &exc)
      {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
      }
    catch (...)
      {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
      }

  return 0;
}
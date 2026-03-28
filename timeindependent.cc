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
#include "equation.h"

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

template <int n>
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
  AssemblyScratchData(const FiniteElement<3> &fe);
  AssemblyScratchData(const AssemblyScratchData &scratch_data);

  FEValues<3> fe_values;

  std::vector<Vector<double>> values_newton;
  std::vector<Vector<double>> values_old;

  std::vector<std::vector<Tensor<1,3>>> gradients_newton;

  std::vector<double> acegen_scratch;
};
 
  struct AssemblyCopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };
  
  void assemble_system();

  void local_assemble_system(
      const typename DoFHandler<3>::active_cell_iterator &cell,
      AssemblyScratchData                                  &scratch,
      AssemblyCopyData                                     &copy_data);
  void copy_local_to_global(const AssemblyCopyData &copy_data);

  void solve();
  double determine_step_length() const;
  void output_results() const;
  void generate_rhs();
  double compute_residual();

  Triangulation<3> triangulation;
  const FESystem<3>    fe;
  DoFHandler<3>    dof_handler;

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

  int max_it;
  double max_multiplier;
  double min_multiplier;
  int optimal_it;
  int newton_iteration;
};

template <int n>
Step3<n>::Step3()
  : fe(FE_Q<3>(1), n)
  , dof_handler(triangulation)
  , n_q_points(QGauss<3>(fe.degree + 1).size())
  , max_it(10)
  , max_multiplier(1.6)
  , min_multiplier(0.5)
  , optimal_it(5)
  , newton_iteration(0)
{}

template <int n>
class RightHandSide : public Function<3>
{
public:
  RightHandSide()
    : Function<3>(n)
  {}

  virtual double value(const Point<3> &p,
                       const unsigned int component = 0) const override;
};


template <int n>
double RightHandSide<n>::value(const Point<3> &p,
                               const unsigned int component) const
{

    double u1 = std::sin(M_PI*p[0]) *
                std::sin(M_PI*p[1]) *
                std::sin(M_PI*p[2]);

    double u2 = std::sin(M_PI*p[0]) *
                std::sin(M_PI*p[1]) *
                std::sin(M_PI*p[2]);

    if (component == 0)
    {
        return 
             3*M_PI*M_PI * u1;
    }
    else if (component == 1)
    {
        return
             6*M_PI*M_PI * u2;
    }

    return 0.0;
}

template <int n>
class ExactSolution : public Function<3>
{
  public:
  ExactSolution()
  : Function<3>(n)
  {}
  
  virtual double value(const Point<3> &p,
                      const unsigned int component = n) const override
  {
    (void)p;
    (void)component;

    double u1 = std::sin(M_PI*p[0])*std::sin(M_PI*p[1])*std::sin(M_PI*p[2]);
    double u2 = std::sin(M_PI*p[0])*std::sin(M_PI*p[1])*std::sin(M_PI*p[2]);

    if (component == 0) {
      return u1;
    }
    else if (component == 1) {
      return u2;
    }

    return 0.0;
  }
};


template <int n>
void Step3<n>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  triangulation.refine_global(
    6
  );

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}

template <int n>
void Step3<n>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  zero_constraints.clear();
  for (int i = 0; i < 6; i++)
  {
    VectorTools::interpolate_boundary_values(dof_handler,
                                                i,
                                                Functions::ZeroFunction<3>(n),
                                                zero_constraints);
  }

  zero_constraints.close();

  nonzero_constraints.clear();

  DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, nonzero_constraints); // x-direction
  DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, nonzero_constraints); // y-direction
  DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, nonzero_constraints); // z-direction

  nonzero_constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints, true);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  oldsolution.reinit(dof_handler.n_dofs());
  newton_iterate.reinit(dof_handler.n_dofs());

  rhs_values.resize(triangulation.n_active_cells(),
                  std::vector<Tensor<1,n>>(n_q_points));
}

template <int n>
Step3<n>::AssemblyScratchData::AssemblyScratchData(const FiniteElement<3> &fe)
  :
  fe_values(fe,
            QGauss<3>(fe.degree + 1),
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
                          std::vector<Tensor<1,3>>(n));

  acegen_scratch.resize(256);
}

template <int n>
Step3<n>::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
    , values_newton(scratch_data.values_newton)
    , values_old(scratch_data.values_old)
    , gradients_newton(scratch_data.gradients_newton)
    , acegen_scratch(scratch_data.acegen_scratch)
  {}

template <int n>
void Step3<n>::assemble_system()
{
  system_matrix = 0;
  system_rhs = 0;

  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &Step3<n>::local_assemble_system,
                  &Step3<n>::copy_local_to_global,
                  AssemblyScratchData(fe),
                  AssemblyCopyData());
}

template <int n>
void Step3<n>::local_assemble_system(
        const typename DoFHandler<3>::active_cell_iterator &cell,
        AssemblyScratchData                                  &scratch_data,
        AssemblyCopyData                                     &copy_data)
{
  const QGauss<3> quadrature_formula(fe.degree + 1);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  //n_q_points    = quadrature_formula.size();

  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);

  copy_data.local_dof_indices.resize(dofs_per_cell);
 
  scratch_data.fe_values.reinit(cell);


 const unsigned int dim=3;

  //Acegen OUTPUT
  //RESIDUAL
  std::vector<Vector<double>> dPsiDu(n_q_points, Vector<double>(n));
  
  std::vector<std::vector<Tensor<1,dim>>> dPsidGradU(
      n_q_points,
      std::vector<Tensor<1,dim>>(n));

  //TANGENT
  std::vector<FullMatrix<double>> dPsiDu2(n_q_points, FullMatrix<double>(n,n));
  
  std::vector<std::vector<std::vector<Tensor<1,dim>>>> dPsidUdGradU(
      n_q_points,
      std::vector<std::vector<Tensor<1,dim>>>(
          n,
          std::vector<Tensor<1,dim>>(n)));

  std::vector<std::vector<std::vector<Tensor<2,dim>>>> dPsidGradU2(
      n_q_points,
      std::vector<std::vector<Tensor<2,dim>>>(
          n,
          std::vector<Tensor<2,dim>>(n)));
  //======= ACEGEN=======

  RightHandSide<n> rhs_function;

  scratch_data.fe_values.reinit(cell);

  copy_data.cell_matrix = 0;
  copy_data.cell_rhs    = 0;

  scratch_data.fe_values.get_function_gradients(solution, scratch_data.gradients_newton);
  scratch_data.fe_values.get_function_values(oldsolution, scratch_data.values_old);
  scratch_data.fe_values.get_function_values(solution, scratch_data.values_newton);

  const unsigned int cell_id = cell->active_cell_index();

  //right_hand_side(fe_values.get_quadrature_points(), rhs_values);

  const FEValuesExtractors::Vector displacements (0);

  for (const unsigned int q_index : scratch_data.fe_values.quadrature_point_indices())
    {
    const auto &x_q = scratch_data.fe_values.quadrature_point(q_index);

    equation(scratch_data.acegen_scratch,
            scratch_data.values_newton[q_index],
        scratch_data.gradients_newton[q_index],
        dPsiDu[q_index],
        dPsidGradU[q_index],
        dPsiDu2[q_index],
        dPsidUdGradU[q_index],
        dPsidGradU2[q_index] // 2, 2, 3, 3
      );

      const auto &sd = scratch_data;

      for (const unsigned int i : scratch_data.fe_values.dof_indices())
      {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        
          for (const unsigned int j : scratch_data.fe_values.dof_indices())
          {
            const unsigned int component_j = fe.system_to_component_index(j).first;

            copy_data.cell_matrix(i, j) +=
                (sd.fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    dPsidGradU2[q_index][component_i][component_j] *           // dPsi/d2(grad u)
                sd.fe_values.shape_grad(j, q_index)    // grad phi_j(x_q)
                + 
                sd.fe_values.shape_value(i, q_index) * // phi_i(x_q)
                    dPsiDu2[q_index][component_i][component_j] *           // dPsi/d2( u)
                sd.fe_values.shape_value(j, q_index)     // phi_j(x_q)
                +
                sd.fe_values.shape_value(i, q_index) * //  phi_i(x_q)
                    dPsidUdGradU[q_index][component_i][component_j] *           // dPsi/d(grad u)du
                sd.fe_values.shape_grad(j, q_index)     // grad phi_j(x_q)
                +
                sd.fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    dPsidUdGradU[q_index][component_i][component_j] *           // dPsi/d(grad u)du
                sd.fe_values.shape_value(j, q_index)     // phi_j(x_q)
                ) *
                sd.fe_values.JxW(q_index);           // dx
        }

                  // number of degrees of freedom is equql to number "geometric dof" = support points
          // function fe.system_to_component_index(i) returns a pair:
          // first is the corresponding component of vector system
          // second is the index of support point
          // shape value known about it through index i
          // ted tim vzdycky zredukuju dimenzi objektu a prevedu to na skalarni pripad

            copy_data.cell_rhs(i) -= (sd.fe_values.shape_value(i, q_index) * 
                            dPsiDu[q_index][component_i]

                        + sd.fe_values.shape_grad(i, q_index) * dPsidGradU[q_index][component_i]
                            - sd.fe_values.shape_value(i,q_index) * rhs_function.value(x_q, component_i)) *
                            sd.fe_values.JxW(q_index);
      }
    }
  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int n>
void Step3<n>::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    zero_constraints.distribute_local_to_global(
      copy_data.cell_matrix,
      copy_data.cell_rhs,
      copy_data.local_dof_indices,
      system_matrix,
      system_rhs);
  }

template <int n>
double Step3<n>::determine_step_length() const
{
  return 1.0;
}

template <int n>
void Step3<n>::solve()
{
  SolverControl            solver_control(5000, 1e-6 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);

  solver.solve(system_matrix, newton_iterate, system_rhs, preconditioner);

  zero_constraints.distribute(newton_iterate);

  const double alpha = determine_step_length();

  solution.add(alpha, newton_iterate);

  zero_constraints.distribute(solution);

  std::cout << solver_control.last_step()
            << " iterations needed to obtain convergence." << std::endl;
}

template <int n>
void Step3<n>::output_results() const
{
  DataOut<3> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  const std::string filename = "output/solution.vtu";
  std::ofstream output(filename);
  data_out.write(output, DataOutBase::vtu);

  std::cout << "Output written to " << filename << std::endl;

  std::ofstream pvd_output ("solution.pvd");
}

template <int n>
double Step3<n>::compute_residual()
{
    ExactSolution<n> exact_solution;

    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                        solution,
                                        exact_solution,
                                        difference_per_cell,
                                        QGauss<3>(fe.degree + 1),
                                        VectorTools::L2_norm);
    const double L2_error =
    VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);
    return L2_error;
}

template <int n>
void Step3<n>::run()
{

  make_grid();
  setup_system();

  std::cout << n_q_points << std::endl;

ExactSolution<n> exact_solution0;

  VectorTools::project(dof_handler,
                      zero_constraints,
                      QGauss<3>(fe.degree + 1),
                      exact_solution0,
                      solution);

  zero_constraints.distribute(solution);

  //std::cout << "The difference is " << L2_error << std::endl;

  output_results();

  newton_iteration = 0;
  double residual_norm = 1.0;
  while (residual_norm > 1e-10 && newton_iteration < max_it)
  {
    //newton_iterate = 0.0;
    assemble_system();
    residual_norm = system_rhs.l2_norm();

    std::cout << "The norm of our solution is: " << residual_norm << std::endl;

    solve();

    newton_iteration++;
  }

    std::cout << "Chyba: " << 100*compute_residual() / solution.linfty_norm() << "%" << std::endl;

    output_results();
}



int main()
{
    using namespace dealii;
    try
      {
        MultithreadInfo::set_thread_limit();
  
        Step3<2> double_ditch;
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
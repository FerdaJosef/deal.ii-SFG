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
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include "nonlinear.h"
#include <cmath>

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

using namespace dealii;

class RightHandSide : public Function<3>
{
public:
  RightHandSide()
    : Function()
    , period(0.2)
  {}

  virtual double value(const Point<3>  &p,
                       const unsigned int component = 0) const override;
  
  private:
    const double period;
};

double RightHandSide::value(const Point<3> &p,
                                 const unsigned int /*component*/) const
{
  const double t = this->get_time();
  const double S = std::pow(std::sin(numbers::PI * p[0]), 2)
                 * std::pow(std::sin(numbers::PI * p[1]), 2)
                 * std::pow(std::sin(numbers::PI * p[2]), 2);
  const double f = S * std::exp(-6*numbers::PI*numbers::PI*t);

  /*const double f = (-1.0 + 12.0 * numbers::PI * numbers::PI) * u_exact
                   + u_exact * u_exact; */

  return f;
}

class BoundaryValues : public Function<3>
{
public:
  virtual double value(const Point<3>  &p,
                       const unsigned int component = 0) const override;
};

double BoundaryValues::value(const Point<3> &p,
                                  const unsigned int /*component*/) const
{
  //const double time = this->get_time();

  return 0*p[0];
  //return 10*std::cos(M_PI*p[0])*std::cos(M_PI*p[1]);
  //return std::cos(p[0] - 11*p[1] + 4*p[2])*std::exp(p[0]*p[1]*p[2]);
}

class Step3
{
public:
  Step3();

  void run();


private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;
  double determine_step_length() const;
  void write_maxima();
  double compute_residual();

  Triangulation<3> triangulation;
  const FE_Q<3>    fe;
  DoFHandler<3>    dof_handler;

  AffineConstraints<double> constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> oldsolution;
  Vector<double> newton_iterate;
  Vector<double> solution;
  Vector<double> system_rhs;

  double time;
  double final_time;
  double delta_t;
  unsigned int timestep_number;

  std::vector<double> time_values;
  std::vector<double> maxima;
};

class ExactSolution : public Function<3>
{
  public:
  ExactSolution(const unsigned int n_components = 1, const double time = 0.)
  : Function<3>(n_components, time)
  {}
  
  virtual double value(const Point<3> &p,
                      const unsigned int /*component*/ = 0) const override
  {
    const double t = this->get_time();

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d{0.5, 0.05};

    double sample;
    sample = d(gen);
    //return 10*std::cos(M_PI*p[0])*std::cos(M_PI*p[1]);
    return sample;

  }
};

class InitialValues : public Function<3>
{
  public:
    InitialValues(const unsigned int n_components = 1, const double time = 0.)
      : Function<3>(n_components, time)
    {}
  
    virtual double value(const Point<3>  &p,
                        const unsigned int component = 0) const override
    {
      return ExactSolution(1, this->get_time()).value(p, component);
    }
};

Step3::Step3()
  : fe(/* polynomial degree = */ 1)
  , dof_handler(triangulation)
  , time(0.0)
  , final_time(1.0)
  , delta_t(0.1)
  , timestep_number(0)
{}

void Step3::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
 
  triangulation.refine_global(4);

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}

void Step3::write_maxima()
{
  std::ofstream out("evolution.txt");
  for (unsigned int i = 0; i < time_values.size(); ++i)
  {
      out << time_values[i] << " " << maxima[i] << "\n";
  }
}

double Step3::compute_residual()
{
  Vector<float> difference_per_cell(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    ExactSolution(1, time),
                                    difference_per_cell,
                                    QGauss<3>(fe.degree + 1),
                                    VectorTools::L2_norm);
  const double L2_error =
  VectorTools::compute_global_error(triangulation,
                                    difference_per_cell,
                                    VectorTools::L2_norm);
  return L2_error;
}

void Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  
  constraints.clear();
  
  //DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints); // x-direction
  //DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints); // y-direction
  //DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, constraints); // z-direction

  VectorTools::interpolate_boundary_values(dof_handler,
                                          types::boundary_id(0),
                                          BoundaryValues(),
                                          constraints);

  constraints.close();

  std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  oldsolution.reinit(dof_handler.n_dofs());
  newton_iterate.reinit(dof_handler.n_dofs());
}



void Step3::assemble_system()
{
  const QGauss<3> quadrature_formula(fe.degree + 1);
  FEValues<3> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);


  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  system_matrix = 0;
  system_rhs = 0;

  RightHandSide rhs_function;
  rhs_function.set_time(this->time);

 const unsigned int dim=3;
  //======= ACEGEN input
  std::vector<double >  values_old(     quadrature_formula.size());
  std::vector<double >  values_newton(quadrature_formula.size());
  std::vector<Tensor<1, dim>>  gradients_newton(quadrature_formula.size());
  std::vector<double> acegen_scratch(256);


  //Acegen OUTPUT
  //RESIDUAL
  std::vector<double> dPsiDU(quadrature_formula.size());
  std::vector<Tensor<1, dim>> dPsiDgradU(quadrature_formula.size());

  //TANGENT
  std::vector<double> dPsiDU2(quadrature_formula.size());
  std::vector<Tensor<1, dim>> dPsiDUdgradu(quadrature_formula.size());
  std::vector<Tensor<2, dim>> dPsiDgradu2(quadrature_formula.size());
  //======= ACEGEN======

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      //fe_values.get_function_gradients(solution, gradients_newton);
      fe_values.get_function_values(oldsolution, values_old);
      fe_values.get_function_values(solution, values_newton);
      fe_values.get_function_gradients(solution, gradients_newton);


      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
        laplace_nonlinear(acegen_scratch,
                &(values_newton[q_index]),
            &(values_old[q_index]),
            gradients_newton[q_index],
            &(dPsiDU[q_index]),
            dPsiDgradU[q_index],
            &dPsiDU2[q_index],
            dPsiDgradu2[q_index],
            dPsiDUdgradu[q_index],
            &this->delta_t
          );

          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    dPsiDgradu2[q_index] *           // dPsi/d(grad u)
                 fe_values.shape_grad(j, q_index)     // grad phi_j(x_q)
                 + 
                 fe_values.shape_value(i, q_index) * // phi_i(x_q)
                    dPsiDU2[q_index] *           // dPsi/d2( u)
                 fe_values.shape_value(j, q_index)     // phi_j(x_q)
                 +
                 fe_values.shape_value(i, q_index) * //  phi_i(x_q)
                    dPsiDUdgradu[q_index] *           // dPsi/d(grad u)du
                 fe_values.shape_grad(j, q_index)     // grad phi_j(x_q)
                 +
                 fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    dPsiDUdgradu[q_index] *           // dPsi/d(grad u)du
                 fe_values.shape_value(j, q_index)     // phi_j(x_q)
                ) *
                 fe_values.JxW(q_index);           // dx

          const auto &x_q = fe_values.quadrature_point(q_index);
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) -= (fe_values.shape_value(i, q_index) * 
                            (dPsiDU[q_index] + rhs_function.value(x_q))
                        + fe_values.shape_grad(i, q_index) * dPsiDgradU[q_index]
                        ) *
                            fe_values.JxW(q_index);
        }
      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                              local_dof_indices,
                                              system_matrix, system_rhs);
    }
  }

double Step3::determine_step_length() const
{
  return 1.0;
}

void Step3::solve()
{
  SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  //solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  solver.solve(system_matrix, newton_iterate, system_rhs, PreconditionIdentity());  

  constraints.distribute(newton_iterate);

  const double alpha = determine_step_length();
  solution.add(alpha, newton_iterate);  

  constraints.distribute(solution);
  //std::cout << solver_control.last_step()
  //          << " CG iterations needed to obtain convergence." << std::endl;
  
  
}

void Step3::output_results() const
{
  static std::vector<std::pair<double, std::string>> times_and_names;

  DataOut<3> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  const std::string filename = "output/solution-" + Utilities::int_to_string(timestep_number) + ".vtu";
  std::ofstream output(filename);
  data_out.write(output, DataOutBase::vtu);

  std::cout << "Output written to " << filename << std::endl;
  times_and_names.push_back(
                {time, filename});

  std::ofstream pvd_output ("solution.pvd");
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}



void Step3::run()
{
  make_grid();
  setup_system();


  VectorTools::project(dof_handler,
                        constraints,
                        QGauss<3>(fe.degree + 1),
                        ExactSolution(1, 0.0),
                        solution);

  constraints.distribute(solution);

  maxima.push_back(solution.linfty_norm());
  time_values.push_back(time);
                
  //std::cout << "The difference is " << L2_error << std::endl;

  output_results();

  std::cout << time << std::endl;
  while (time < final_time)
  {
    time+=delta_t;
    timestep_number++;

    oldsolution = solution;
    double residual_norm = 1.0;
    unsigned int newton_iteration = 0;
    while (residual_norm > 1e-10 && newton_iteration < 200)
    {
      newton_iterate = 0.0;
      assemble_system();
      residual_norm = system_rhs.l2_norm();

      std::cout << "The norm of our solution is: " << residual_norm << std::endl;

      solve();

      newton_iteration++;
    }
    maxima.push_back(solution.linfty_norm());
    time_values.push_back(time);

    std::cout << "The maximum is: " << maxima.back() << std::endl;
    output_results();
  }
  write_maxima();
}



int main()
{
  Step3 laplace_problem;
  laplace_problem.run();
  return 0;
}

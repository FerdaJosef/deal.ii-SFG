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
#include "laplace.h"
#include <filesystem>

using namespace dealii;

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
 
      FEValues<3>     fe_values;
 
      std::vector<double> rhs_values;
      BoundaryValues<3> boundary_values;
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
  void output_results(/*const unsigned int cycle*/) const;

  Triangulation<3> triangulation;
  const FE_Q<3>    fe;
  DoFHandler<3>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> oldsolution;
  Vector<double> newton_iterate;
  Vector<double> solution;
  Vector<double> system_rhs;

  double       time;
  double       time_step;
  unsigned int timestep_number;

  unsigned int iterations;
};

class BoundaryValues : public Function<3>
{
public:
  virtual double value(const Point<3>  &p,
                       const unsigned int component = 0) const override;
  BoundaryValues() : Function<3>() {}
};

double BoundaryValues::value(const Point<3> &p,
                                  const unsigned int /*component*/) const
{
  const double t = this->get_time();

  const double tol = 1e-10; // tolerance
  /*
  if (std::fabs(p[0] + 1.0) < tol)
    return 0.0;

  if (std::fabs(p[0] - 1.0) < tol)
    return 1.0;

  if (std::fabs(p[1] + 1.0) < tol)
    return 0.0;

  if (std::fabs(p[1] - 1.0) < tol)
    return 1.0;

  if (std::fabs(p[2] + 1.0) < tol)
    return 0.0;

  if (std::fabs(p[2] - 1.0) < tol)
    return 1.0;
  */

  return std::cos(5.0*M_PI*p[0])*std::exp(-2500.0*(t - 0.05)*(t - 0.05));
}

Step3::Step3()
  : fe(/* polynomial degree = */ 1)
  , dof_handler(triangulation)
{}



void Step3::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(4);

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}

void Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  oldsolution.reinit(dof_handler.n_dofs());
  newton_iterate.reinit(dof_handler.n_dofs());
}

void Step3::assemble_system()
{
  WorkStream::run(dof_handler.begin_active(),
                  dof_handler.end(),
                  *this,
                  &Step3::local_assemble_system,
                  &Step3::copy_local_to_global,
                  AssemblyScratchData(fe),
                  AssemblyCopyData());
}

Step3::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<3> &fe)
    : fe_values(fe,
                QGauss<3>(fe.degree + 1),
                update_values | update_gradients | update_quadrature_points |
                  update_JxW_values)
  {}

Step3::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  {}

void Step3::assemble_local_system(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AssemblyScratchData                                  &scratch_data,
    AssemblyCopyData                                     &copy_data)
{
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  // I dont think we need this! We already have it in AssemblyCopyData

  //FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  //Vector<double>     cell_rhs(dofs_per_cell);

  // But we do have to reinit after every step, right?

  //std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  // I left off here!!!!!!!

  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);
  copy_data.local_dof_indices.resize(dofs_per_cell);
  scratch_data.fe_values.reinit(cell);


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
  //======= ACEGEN=======

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.get_function_gradients(newton_iterate, gradients_newton);
      fe_values.get_function_values(oldsolution, values_old);
      fe_values.get_function_values(newton_iterate, values_newton);


      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
        laplace(acegen_scratch,
                &(values_newton[q_index]),
            &(values_old[q_index]),
            gradients_newton[q_index],
            &(dPsiDU[q_index]),
            dPsiDgradU[q_index],
            &dPsiDU2[q_index],
            dPsiDgradu2[q_index],
            dPsiDUdgradu[q_index],
            &time_step, &time
          );

          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              copy_data.cell_matrix(i, j) +=
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

          // Toto je rovnice pro konečné prvky        
          for (const unsigned int i : fe_values.dof_indices())
            copy_data.cell_rhs(i) += (fe_values.shape_value(i, q_index) * 
                            dPsiDU[q_index]
                        + fe_values.shape_grad(i, q_index) * dPsiDgradU[q_index]
                        ) *
                            fe_values.JxW(q_index);
        }
      cell->get_dof_indices(copy_data.local_dof_indices);

      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(copy_data.local_dof_indices[i],
                            copy_data.local_dof_indices[j],
                            copy_data.cell_matrix(i, j));

      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  std::map<types::global_dof_index, double> boundary_values;

  BoundaryValues BV;
  BV.set_time(time);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           types::boundary_id(0),
                                           BV,
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}



void Step3::solve()
{
  SolverControl            solver_control(1000, 1e-06*system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  iterations = solver_control.last_step();
  std::cout << iterations
            << " CG iterations needed to obtain convergence." << " " << system_rhs.l2_norm() << std::endl;
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
  DataOutBase::write_pvd_record (pvd_output, times_and_names);
}


/*
void Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}
*/

void Step3::run()
{
  std::filesystem::create_directory("output");

  time = 0.0;
  time_step = 0.001;
  timestep_number = 0;
  const double final_time = 0.1; // for example

  make_grid();
  setup_system();

  oldsolution = 0.0;
  solution = 0.0;

  double max_multiplier = 1.6;
  double min_multiplier = 0.5;
  double optimal_it = 30;
  double dt_max = final_time / 10;
  
  while (time < final_time)
  {
      time += time_step;
      ++timestep_number;
      std::cout << "Time step " << timestep_number << " at t=" << time << std::endl;

      assemble_system();
      solve();

      if (iterations <= optimal_it) {
        // self.dt.assign(min(dt_max, float(self.dt) * min(1.6, optimal_it / (iterations + 0.001))));
        time_step = std::min(dt_max, time_step*std::min(max_multiplier, optimal_it/(iterations + 0.001)));
      }
      else {
        // self.dt.assign(float(self.dt) * max(0.5, optimal_it / iterations));
        time_step = time_step*std::max(min_multiplier, optimal_it / iterations);
      }

      if (time + time_step > final_time) {
        time_step = final_time - time;
      }
      
      if (time_step < 1e-6) {
        break;
      }

      output_results();

      oldsolution = solution;
  }
}


int main()
{
  Step3 laplace_problem;
  laplace_problem.run();

  return 0;
}

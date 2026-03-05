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

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include "equation.h"

using namespace dealii;

class BoundaryValues : public Function<3>
{
public:
  BoundaryValues();
  virtual double value(const Point<3>  &p,
                       const unsigned int component = 0) const override;
};

BoundaryValues::BoundaryValues()
  : Function<3>(2)
{}

double BoundaryValues::value(const Point<3> &p,
                                  const unsigned int /*component*/) const
{
  return p[1];
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

  Triangulation<3> triangulation;
  const FESystem<3>    fe;
  DoFHandler<3>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> oldsolution;
  Vector<double> newton_iterate;
  Vector<double> solution;
  Vector<double> system_rhs;
};


Step3::Step3()
  : fe(FE_Q<3>(1), 2)
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
  const QGauss<3> quadrature_formula(fe.degree + 1);
  FEValues<3> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


 const unsigned int dim=3;
  //======= ACEGEN input
  std::vector<Vector<double>> values_newton(n_q_points, Vector<double>(2));

  std::vector<Vector<double>> values_old(n_q_points, Vector<double>(2));

  std::vector<std::vector<Tensor<1, dim>>> gradients_newton(n_q_points, std::vector<Tensor<1, dim>>(2));
// USING 
  std::vector<double> acegen_scratch(256);


  //Acegen OUTPUT
  //RESIDUAL
  std::vector<Vector<double>> dPsiDu(n_q_points, Vector<double>(2));
  
  std::vector<std::vector<Tensor<1,dim>>> dPsidGradU(
      n_q_points,
      std::vector<Tensor<1,dim>>(2));

  //TANGENT
  std::vector<FullMatrix<double>> dPsiDu2(n_q_points, FullMatrix<double>(2,2));
  
  std::vector<std::vector<std::vector<Tensor<1,dim>>>> dPsidUdGradU(
      n_q_points,
      std::vector<std::vector<Tensor<1,dim>>>(
          2,
          std::vector<Tensor<1,dim>>(2)));

  std::vector<std::vector<std::vector<Tensor<2,dim>>>> dPsidGradU2(
      n_q_points,
      std::vector<std::vector<Tensor<2,dim>>>(
          2,
          std::vector<Tensor<2,dim>>(2)));
  //======= ACEGEN=======

  double delta_t=0.01;


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      fe_values.get_function_gradients(solution, gradients_newton);
      fe_values.get_function_values(oldsolution, values_old);
      fe_values.get_function_values(solution, values_newton);


      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
        equation(acegen_scratch,
                values_newton[q_index],
            gradients_newton[q_index],
            dPsiDu[q_index],
            dPsidGradU[q_index],
            dPsiDu2[q_index],
            dPsidUdGradU[q_index],
            dPsidGradU2[q_index] // 2, 2, 3, 3 
          );

          for (const unsigned int i : fe_values.dof_indices())
            {
            const unsigned int component_i = fe.system_to_component_index(i).first;
            
            // number of degrees of freedom is equql to number "geometric dof" = support points
            // function fe.system_to_component_index(i) returns a pair:
            // first is the corresponding component of vector system
            // second is the index of support point
            // shape value known about it through index i
            // ted tim vzdycky zredukuju dimenzi objektu a prevedu to na skalarni pripad

              cell_rhs(i) += (fe_values.shape_value(i, q_index) * 
                              dPsiDu[q_index][component_i]

                          + fe_values.shape_grad(i, q_index) * dPsidGradU[q_index][component_i]
                          ) *
                              fe_values.JxW(q_index);
              }
          for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i = fe.system_to_component_index(i).first;
            
              for (const unsigned int j : fe_values.dof_indices())
              {
                const unsigned int component_j = fe.system_to_component_index(j).first;

                cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                        dPsidGradU2[q_index][component_i][component_j] *           // dPsi/d2(grad u)
                    fe_values.shape_grad(j, q_index)    // grad phi_j(x_q)
                    + 
                    fe_values.shape_value(i, q_index) * // phi_i(x_q)
                        dPsiDu2[q_index][component_i][component_j] *           // dPsi/d2( u)
                    fe_values.shape_value(j, q_index)     // phi_j(x_q)
                    +
                    fe_values.shape_value(i, q_index) * //  phi_i(x_q)
                        dPsidUdGradU[q_index][component_i][component_j] *           // dPsi/d(grad u)du
                    fe_values.shape_grad(j, q_index)     // grad phi_j(x_q)
                    +
                    fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                        dPsidUdGradU[q_index][component_i][component_j] *           // dPsi/d(grad u)du
                    fe_values.shape_value(j, q_index)     // phi_j(x_q)
                    ) *
                    fe_values.JxW(q_index);           // dx
            }
          }
        }
      cell->get_dof_indices(local_dof_indices);

      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           types::boundary_id(0),
                                           BoundaryValues(),
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}



void Step3::solve()
{
  SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  std::cout << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}



void Step3::output_results() const
{
  DataOut<3> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  const std::string filename = "solution.vtk";
  std::ofstream     output(filename);
  data_out.write_vtk(output);
  std::cout << "Output written to " << filename << std::endl;
}



void Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}



int main()
{
  Step3 laplace_problem;
  laplace_problem.run();

  return 0;
}

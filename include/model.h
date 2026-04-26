#ifndef DOUBLE_DITCH
#define DOUBLE_DITCH

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
 
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

template <int dim, int n>
class Step3
{
public:
  Step3(ParameterHandler &);

  void run();


private:

  ParameterHandler &prm;

  void make_grid();
  void setup_system();
  void parse_parameters();

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
  bool time_step_update();
  double determine_step_length() const;
  void output_results() const;
  void generate_rhs();
  double compute_residual();
  void make_timestep();

  Triangulation<dim> triangulation;
  const FESystem<dim>    fe;
  DoFHandler<dim>    dof_handler;

  AffineConstraints<double> constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> oldsolution;
  Vector<double> newton_iterate;
  Vector<double> solution;
  Vector<double> system_rhs;

  const unsigned int n_q_points;
  std::vector<std::vector<Tensor<1,n>>> rhs_values;

  unsigned int n_refinements;
  double left_lim;
  double right_lim;

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

#endif //DOUBLE_DITCH
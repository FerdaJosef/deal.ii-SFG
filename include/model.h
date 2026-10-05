#ifndef DOUBLE_DITCH
#define DOUBLE_DITCH

// --- Mesh & Geometry ---
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

// --- Degrees of Freedom (DoFs) ---
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

// --- Finite Elements & Quadrature ---
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/base/quadrature_lib.h>

// --- Linear Algebra ---
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>

// --- Linear Solvers & Preconditioners ---
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>

// --- Parallelization (MPI) ---
#include <deal.II/base/multithread_info.h>

// --- Numerics & Tools ---
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>

// --- Standard Library ---
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <filesystem>

// --- External Headers ---
#include "equation.h"
#include "RandomField.h"

// --- MPI & PETSc ---
#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

// --- Adaptive Mesh Refinement

using namespace dealii;

template <int dim, int n>
class Step3
{
public:
  Step3(ParameterHandler &param);
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  bool time_step_update();
  double determine_step_length() const;
  void output_results() const;
  void make_timestep();

  // --- MPI Controls ---
  MPI_Comm mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
  ConditionalOStream pcout;

  ParameterHandler &prm;
  mutable TimerOutput computing_timer;

  // --- Physics & Random Field ---
  RandomField<dim, n> random_field;

  // --- Mesh & FE Data ---
  Triangulation<dim> triangulation;
  const FESystem<dim> fe;
  DoFHandler<dim> dof_handler;
  AffineConstraints<double> constraints;

  // --- PETSc Linear Algebra ---
// --- PETSc Linear Algebra ---
  PETScWrappers::MPI::SparseMatrix system_matrix;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  // Non-ghosted vectors (for math, projections, and linear algebra)
  PETScWrappers::MPI::Vector distributed_solution;
  PETScWrappers::MPI::Vector distributed_old_solution;
  PETScWrappers::MPI::Vector newton_iterate;
  PETScWrappers::MPI::Vector system_rhs;

  // Ghosted vectors (read-only for assembly evaluation and output)
  PETScWrappers::MPI::Vector solution;
  PETScWrappers::MPI::Vector oldsolution;

  const unsigned int n_q_points;

  // --- Simulation Parameters ---
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
  double dt_max;
  double dt_min;
  int newton_iteration;
  int solver_iteration;

  double linear_residual;
  int max_linear_iteration;
};

#endif // DOUBLE_DITCH
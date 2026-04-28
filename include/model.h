#ifndef DOUBLE_DITCH
#define DOUBLE_DITCH

// --- Mesh & Geometry ---
#include <deal.II/grid/tria.h>                // Basic triangulation classes
#include <deal.II/grid/grid_generator.h>      // Standard shapes (hyper_cube, etc.)
#include <deal.II/grid/grid_tools.h>          // Transformations, periodicity, finding cells

// --- Degrees of Freedom (DoFs) ---
#include <deal.II/dofs/dof_handler.h>         // Manages DoF distribution on the mesh
#include <deal.II/dofs/dof_tools.h>           // High-level DoF operations (constraints, etc.)
#include <deal.II/dofs/dof_renumbering.h>     // Reordering DoFs (Cuthill-McKee, etc.) for speed

// --- Finite Elements & Quadrature ---
#include <deal.II/fe/fe_q.h>                  // Lagrange finite elements
#include <deal.II/fe/fe_system.h>             // Composing vector-valued elements (n-variables)
#include <deal.II/fe/fe_values.h>             // Shape function evaluation at quadrature points
#include <deal.II/fe/mapping_q.h>             // Maps reference cell to real cell geometry
#include <deal.II/base/quadrature_lib.h>      // Standard Gauss-Legendre quadrature rules

// --- Linear Algebra (LAC) ---
#include <deal.II/lac/vector.h>               // Simple vector class
#include <deal.II/lac/full_matrix.h>          // Dense matrices (for local assembly)
#include <deal.II/lac/sparse_matrix.h>        // Compressed Row Storage (for global system)
#include <deal.II/lac/dynamic_sparsity_pattern.h> // Building the sparse matrix structure
#include <deal.II/lac/affine_constraints.h>   // Handling Dirichlet, Periodic, and Hanging node constraints

// --- Linear Solvers & Preconditioners ---
#include <deal.II/lac/solver_gmres.h>         // GMRES iterative solver
#include <deal.II/lac/solver_cg.h>            // Conjugate Gradient solver
#include <deal.II/lac/precondition.h>         // Basic preconditioners (Jacobi, etc.)

// --- Parallelization & Multi-threading ---
#include <deal.II/base/work_stream.h>         // Thread-safe assembly management
#include <deal.II/base/multithread_info.h>    // Controls the number of threads used

// --- Numerics & Tools ---
#include <deal.II/numerics/vector_tools.h>    // Interpolation, boundary conditions
#include <deal.II/numerics/matrix_tools.h>    // Assembly of standard matrices (Mass, Laplace)
#include <deal.II/numerics/data_out.h>        // Writing .vtu/.vtk files for Paraview
#include <deal.II/numerics/error_estimator.h> // Kelly error estimator (if needed for refinement)
#include <deal.II/base/parameter_handler.h>   // Reading inputs from .prm files
#include <deal.II/base/conditional_ostream.h> // Prevents output spam in parallel runs
#include <deal.II/base/timer.h>               // Timer

// --- Standard Library & External ---
#include <fstream>                            // File streams for output
#include <iostream>                           // Console output
#include <math.h>                             // Basic math functions
#include <random>                             // C++ random number generators
#include "equation.h"                         // Your AceGen-generated physics
#include "RandomField.h"                      // Your custom Monte Carlo field

using namespace dealii;
template <int dim, int n>
class Step3
{
public:
  Step3(ParameterHandler &);

  void run();


private:



  ParameterHandler &prm;

  mutable TimerOutput computing_timer;

  RandomField<dim, n> random_field;

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

  double linear_residual;

  std::vector<unsigned int> component_indices;
};

#endif //DOUBLE_DITCH

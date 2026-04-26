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

// Parameter handling

#include <deal.II/base/parameter_handler.h>

#include <random>
#include <filesystem>

//#include "AceGen/2D_3omega/equation_2_3.h"

#include "InitialValues.h"
#include "model.h"
#include "parameters.h"


int main()
{
    using namespace dealii;

    try
      {
        MultithreadInfo::set_thread_limit();

        ParameterHandler prm;
        ParameterReader  param(prm);


        param.read_parameters("double_ditch.prm");

        prm.enter_subsection("Output parameters");
          std::string out_file = prm.get("Output filename"); 
        prm.leave_subsection();

        std::filesystem::path p(out_file);
        if (p.has_parent_path()) 
        {
            std::filesystem::create_directories(p.parent_path());
        }
  
        Step3<3, 2> double_ditch(prm);
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
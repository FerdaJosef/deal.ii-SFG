#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>

#include <filesystem>
#include <iostream>

#include "model.h"
#include "parameters.h"

int main(int argc, char *argv[])
{
  using namespace dealii;

  try
    {
      // 1. Initialize MPI and PETSc. Setting the 3rd argument to 1 limits 
      // each MPI process to 1 thread since WorkStream was removed.
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      ParameterHandler prm;
      ParameterReader  param(prm);

      param.read_parameters("double_ditch.prm");

      prm.enter_subsection("Output parameters");
      std::string out_file = prm.get("Output filename"); 
      prm.leave_subsection();

      // 2. Restrict directory creation to MPI Rank 0 to avoid filesystem collisions
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::filesystem::path p(out_file);
          if (p.has_parent_path()) 
            {
              std::filesystem::create_directories(p.parent_path());
            }
        }

      Step3<2, 4> double_ditch(prm);
      double_ditch.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
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
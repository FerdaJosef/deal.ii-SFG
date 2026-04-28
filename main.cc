
// Parameter handling

#include <deal.II/base/parameter_handler.h>

#include <filesystem>

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
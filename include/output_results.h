#include "model.h"

template <int dim, int n>
void Step3<dim, n>::output_results() const
{
  static std::vector<std::pair<double, std::string>> times_and_names;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  prm.enter_subsection("Output parameters");
 
  const std::string output_filename = prm.get("Output filename");
  data_out.parse_parameters(prm);
 
  prm.leave_subsection();

  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  const std::string filename = output_filename + Utilities::int_to_string(timestep_number) + ".vtu";
  std::ofstream output(filename);
  data_out.write(output, DataOutBase::vtu);

  std::cout << "Output written to " << filename << std::endl;
  times_and_names.push_back(
                {time, filename});

  std::ofstream pvd_output (output_filename);
  DataOutBase::write_pvd_record(pvd_output, times_and_names);
}
  #include "model.h"
  
  class ParameterReader : public EnableObserverPointer
  {
  public:
    ParameterReader(ParameterHandler &);
    void read_parameters(const std::string &);
 
  private:
    void              declare_parameters();
    ParameterHandler &prm;
  };
 
  ParameterReader::ParameterReader(ParameterHandler &paramhandler)
    : prm(paramhandler)
  {}
 
 
  void ParameterReader::declare_parameters()
  {
    prm.enter_subsection("General information");
    {
      prm.declare_entry("Dimension", "2", Patterns::Integer(0), "Physical dimension of our simulation");
      prm.declare_entry("Number of variables", "2", Patterns::Integer(0), "Number of variables in our problem");
    }
    prm.leave_subsection();


    prm.enter_subsection("Mesh & geometry parameters");
    {
        prm.declare_entry("Number of refinements",
                        "6",
                        Patterns::Integer(0),
                        "Number of global mesh refinement steps "
                        "applied to initial coarse grid");

        prm.declare_entry("Mesh left limit",
                        "-25",
                        Patterns::Double(),
                        "Our mesh will have [-left, right] size");

        prm.declare_entry("Mesh right limit",
                        "25",
                        Patterns::Double(),
                        "Our mesh will have [-left, right] size");
    }
    prm.leave_subsection();

    prm.enter_subsection("Temporal parameters");

      prm.declare_entry("Initial time",
                      "6",
                      Patterns::Double(),
                      "Initial time at the beginning of our simulation");

      prm.declare_entry("Final time",
                      "500",
                      Patterns::Double(),
                      "Final time of our simulation");

      prm.declare_entry("Initial time step",
                      "1e-5",
                      Patterns::Double(0),
                      "Initial time step at the beginning of our presentation");

      prm.declare_entry("Max time step",
                        "10.0",
                        Patterns::Double(0),
                        "Maximal possible time step");

      prm.declare_entry("Min time step",
                        "1e-7",
                        Patterns::Double(0),
                        "Minimal possible time step");

      prm.declare_entry("Time step number",
                        "0",
                        Patterns::Integer(0),
                        "Initial time step number (time step 1, time step 2, etc)");

      prm.declare_entry("Max iterations",
                        "10",
                        Patterns::Integer(0),
                        "Maximal possible number of Newton iterations");

      prm.declare_entry("Optimal iterations",
                        "6",
                        Patterns::Integer(0),
                        "Optimal number of Newton iterations");

      prm.declare_entry("Max multiplier",
                        "2.0",
                        Patterns::Double(0),
                        "Maximal multiplier of our time step");

      prm.declare_entry("Min multiplier",
                        "0.5",
                        Patterns::Double(0),
                        "Minimal multiplier of our time step");

    prm.leave_subsection();

    prm.enter_subsection("Output parameters");
    {
      prm.declare_entry("Output filename",
                        "solution",
                        Patterns::Anything(),
                        "Name of the output file (without extension)");
 
      DataOutInterface<1>::declare_parameters(prm);
    }
    prm.leave_subsection();
  }
 
 
  void ParameterReader::read_parameters(const std::string &parameter_file)
  {
    declare_parameters();
 
    prm.parse_input(parameter_file);
  }
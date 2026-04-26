deal.II-SFG Simulation Framework
Technical Documentation & Implementation Guide

This framework provides a structured pipeline for running stochastic phase-field simulations using deal.II and AceGen.
🚀 Setup & Execution Pipeline
Step I: Repository Initialization

Clone the repository and enter the project root:
Bash

git clone <repository-url>
cd deal.ii-SFG

Step II: Model Customization

The model logic is decentralized into specific source files to optimize compilation times.

    Initial Conditions: Modify source/Model/InitialValues.cc.

    Right-Hand Side (RHS): Define source terms and noise in source/Model/RandomField.cc.

    Templates: You MUST append the explicit template instantiation at the end of each .cc file (InitialValues.cc, RandomField.cc, and model.cc):
    C++

    template class RandomField<1, 2>;

Step III: Symbolic Differentiation (AceGen)

The framework utilizes AceGen for automated residual and Jacobian generation.

    Place AceGen C-output in source/Model/AceGen/.

    Run the bridge script to format the code for deal.II:
    Bash

    python3 script.py source/Model/AceGen/output.c

Step IV: Configuration

    main.cc: Include the correct model header and ensure the solver is instantiated with the correct templates.

    CMakeLists.txt: Verify that the TARGET_SRC list includes all relevant .cc files.

Step V: Build & Run

Navigate to the build directory (create it if it doesn't exist) and compile:
Bash

mkdir -p build && cd build
cmake ..
make -j$(nproc)
make run

💡 Technical Notes

    Compilation Optimization: The separation of RandomField.cc and InitialValues.cc ensures that changes to the core FEM algorithm in model.cc do not require re-generating the stochastic noise logic or initial state, significantly reducing iterative development time.

    Troubleshooting: If the linker throws an undefined reference error, double-check that the template class line at the end of your .cc files matches the dim and n used in your simulation.

    Execution: If the program hangs at 100% CPU without output, verify the solver convergence tolerances and the scaling of the random field amplitude.

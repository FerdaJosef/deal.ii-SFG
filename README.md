# Project Setup and Usage Guide

This guide explains how to set up, configure, and run the project with a custom model.

---

## 1. Clone the repository

---

## 2. Create / Configure your model

We put models into the folder "source". We already have three models in this folder:
* Model1D_2vars - 1 dimensional code with 2 omega variants
* Model2D_3vars - 2 dimensional code with 3 omega variants,
* Model3D_2vars - 3 dimensional with 2 omega variants.

Our program has three main classes: Step3 (this includes assembly, timestepping, solving, ... ). InitialValues (self-explanatory) and RandomField (for generating a random right hand side). 
Each model folder includes the following:
* InitialValues.cc
* RandomField.cc
* model.cc
* double_ditch.prm
* AceGen folder

Some of these models are explored in the article "Multiwell phase-field model for arbitrarily strong total-spreading case" by J. Kozlík, K. Tůma, et. al.

If we want to create a custom model, we can do:
```
mkdir source/CustomModel
```

We may copy contents of a previous model into our new one and proceed to make changes from there:
```
cp -a source/Model3D_2vars/. source/CustomModel
```

* We will mostly want to change the following:

  * Right-hand side (RHS)
  * Initial conditions
  * Parameters of our problem (initial time step, final time, number of global refinements, ...).
    We can find parameters in the double_ditch.prm file. The structure of this file can be changed, but corresponding
    changes have to be made to "include/parameters.h".

* Ensure that each `.cc` file containing template code ends with the correct explicit instantiation, for example if we want a 2D problem with 2 omega variants:

```cpp
template class Step3<2,2>;
```

Adjust `<dim, n>` according to your model.

---

## 3. Add AceGen output

1. Place your AceGen-generated file into:

```
source/Model/.../AceGen/
```

2. Run the postprocessing script:

```bash
python3 script.py source/Model/.../AceGen/output.c
```

This will:

* remove incompatible includes (`sms.h`)
* replace unsupported syntax (e.g. `Power`)
* generate a usable `equation.h`

---

## 4. Update `main.cc`

Modify `main.cc` to use the correct model and template parameters:

```cpp
Step3<dim, n> simulation;
simulation.run();
```

Make sure `<dim, n>` matches:

* spatial dimension
* number of variables in your model

---

## 5. Update `CMakeLists.txt`

Ensure the correct include paths are set, especially for:

* your model directory
* AceGen-generated headers

Example:

```cmake
set(MODEL_DIR "${CMAKE_SOURCE_DIR}/source/CustomModel")
```

---

## 6. Build the project

Create and enter a build directory:

```bash
mkdir build
cd build
```

Run CMake:

```bash
cmake ..
```

If this does not work, we have to specify a path to where we installed deal.II.

```bash
cmake -DDEAL_II_DIR=/path/to/deal.II ..
```

Then compile:

```bash
make
```

---

## 7. Run the program

```bash
make run
```

---

## Notes

* Output files are written to the `results/` directory (you can choose the name of the folder in double_ditch.prm).
* If compilation fails due to missing C++17 features, ensure your compiler and CMake configuration use:

```cmake
set(CMAKE_CXX_STANDARD 17)
```

* When switching models, you typically need to:

  * update the AceGen file
  * rerun the Python script
  * rebuild the project

---

## Summary Workflow

```text
Clone → Configure Model → Add AceGen → Run Script → Update main/CMake → Build → Run
```

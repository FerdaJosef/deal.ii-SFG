# Project Setup and Usage Guide

This guide explains how to set up, configure, and run the project with a custom model.

---

## 1. Clone the repository

```bash
git clone <your-repo-url>
cd <your-repo>
```

---

## 2. Create / Configure your model

Inside the `source/Model/...` directory:

* Define or modify:

  * **Right-hand side (RHS)**
  * **Initial conditions**

* Ensure that each `.cc` file containing template code ends with the correct explicit instantiation, for example:

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
include_directories(source/Model/YourModel)
```

If deal.II is not found automatically, you may need:

```bash
cmake -DDEAL_II_DIR=/path/to/deal.II ..
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

* Output files are written to the `results/` directory (created automatically if needed).
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

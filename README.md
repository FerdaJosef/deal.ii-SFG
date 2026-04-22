# deal.ii-SFG

Right now, there are 3 working test examples (one dimensional code with 2 omegas in main_1D.cc, two dimensional with 3 omegas in main_2D.cc, three dimensional with 2 omegas in main_3D.cc).
If you want to run main_1D.cc, go into CMakeLists.txt and use: set(TARGET "main_1D")

To run this code, create a folder where you want your output.
mkdir builder
cd builder
cmake ..    (if this fails, you need to provide full path to installed deal.II: 
            cmake -DDEAL_II_DIR=/path/to/installed/deal.II ..)

make run
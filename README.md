# Lid_Driven_Cavity

Modular Lid Driven cavity code written in C++ using the Finite Difference Method. Preprocessing involves providing Initial and Boundary conditions, defining the domain and selecting the solver control schemes. Postprocessing is handled through MATLAB which creates .mp4 videos of the physical parameters for visualization. 

The C++ code at first reads the input data and displays it, then it shows the simulation and convergence, and finally the MATLAB script to create the .mp4 file is executed from within the C++ code itself. 

Results of 2 test cases are provided - one with Reynolds Number(Re) = 100 and other with Re = 1000.

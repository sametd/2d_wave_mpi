# **2D Wave Equation Solution with MPI**
**This code is written as practice. It includes:**

 - Easily modifiable for any equation
 - Derived Data Types 
 - Ghost point exchanges 
 - Output using posix api
 - Output using NetCDF
 - Very robust Gauss-Seidel solver

**Important Note:** This code designed to run with 1,4 or 9 CPU cores.  It scales well with the increasing processor count, for more core count, grid size and processor size should be compatible. Scaling is mostly IO bound.

Example Visualization for the code output:

<img src="example_output.gif?raw=true" width="200px">


#  **2D Wave Equation Solution with MPI**  
**This code is written as practice. It includes:**  
  
- Easily modifiable for any equation  
- Derived Data Types  
- Cartesian Virtual Topology 
- Ghost point exchanges  
- Output using posix api  
- Output using NetCDF  
- Very robust Gauss-Seidel solver  
  
**Important Note:** This code designed to run with 1,4 or 9 CPU cores. It scales well with the increasing processor count, for more core count, grid size and processor size should be compatible. Scaling is mostly IO bound.  
  
Example Visualization for the code output:  
  
<img  src="example_output.gif?raw=true"  width="200px">  
  
**Compilation Instructions:**  
  
mpicc -O3 -march=native waveEq.c -o wave.x -lnetcdf -lm  
  
You can disable NetCDF support by defining write\_netcdf from config.h and then you can comment out include line for the netcdf.h in the waveEq.c. 

Then, you can compile the code with the following command:
  
mpicc -O3 -march=native waveEq.c -o wave.x -lm

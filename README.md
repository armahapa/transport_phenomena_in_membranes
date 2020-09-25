# transport_phenomena_in_membranes
This repository includes a fortran90 code for a coupled problem of diffusion of curvature elastic protein in cellular membrane with the flow of lipid in plane (2D syrface: code name diff_flow_bend.f90) and also a 1-dimensional code for coupling of diffusion and bending (code name: diff_bend_1D.f90) and the readme file to operate the code 
# Code details: diff_flow_bend.f90
Title: Numerical simulation of a square lipid bilayer with curvature induced-proteins
Simulation code by: Arijit Mahapatra, David Saintillan, Padmini Ranganami (UCSD)
For details of the model, refer to:
"Transport phenomena in fluid films with curvature elasticity" by Mahapatra et al.
 (https://arxiv.org/abs/2001.07539)
# Operation procedure
1. Giving parameters value
You can give parameter values in two ways
(a) Give all the parameters in line 24 or 
    -For this case please activate all the derived non dimensional numbers in either 
    line 49 to line 53 or line 56 to 60
(b) Give the values for all the basic non dimensional numbers in line 40
    -For this case please deactivate the parameter value (line 25-34) and deactivate the 
    expression of derived nondimensional number in lines 49-53. Activate the line 56-60

2. Specifying Grid size
Specify grid numbers in line 63 to 65


3. Specify the initial condition for protein patch
	(a) for 1 patch use lct=1.78d-1 activate xa=0, ya=0 and first line of sig expression
	in line 208
	(b)for 2 patches lct=1.26d-1 activate xa=ya=-0.25, xb=yb=0.25 and first 2 lines of 
	sig expression in line 208
	(a)for 4 patches lct=0.89d-1 activate xa=ya=-0.25, xb=yb=0.25, xc=-yc=-0.25, and
	 xb=-yb=0.25 and first 4 lines of sig expression in line 208

4. Data reading 
   For all the parameters for a particular time step (k) data in a particular row indicates 
   data along y direction and data in particular column represents data along x axis. 
   After N (value of grid size along x direction) rows the data repeats for time step
   k+1.
   
# Code details: diff_bend_1D.f90

# Operation procedure
1. Giving parameters value
Give all the parameters in line 30 to line 41
2. Specifying Grid size
Specify grid numbers in line 55 to 56
3. Specify the initial condition for protein patch in line 75
4. For all the parameters for a particular time step (k) data in a particular column indicates 
   data along the length and data in particular row represents data along with increasing time step. 
   
    


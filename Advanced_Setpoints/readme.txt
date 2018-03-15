***************************************************************
*  Preliminary release of Feedforward Motion Control Toolbox  *
***************************************************************
*              for MatLab Release 12 and up                   *
***************************************************************

Paul Lambrechts, June 15, 2004
Eindhoven University of Technology
Faculty of Mechanical Engineering, Control Systems Technology Group
P.O. Box 513, 5600 MB Eindhoven, The Netherlands
Telephone: +31 40 2472839
Fax:       +31 40 2451418
Email: P.F.Lambrechts@tue.nl
***************************************************************

Properties:
This toolbox allows planning of 2nd, 3d and 4th order trajectories 
for single axis motion systems and using them in feedforward 
control.
The planners are limited to symmetrical trajectories for 
point-to-point motions.
Discrete time implementation is fully supported, dealing with 
quantization effects is indicated.
Use of planners and feedforward controllers in Simulink is 
supported.
***************************************************************

Main Contents:
 make4.m      Matlab function to calculate timing for symmetrical 
              4th order profiles (analytical solution)
 profile4.m   Matlab function to calculate symmetrical 4th order 
              profiles from timing info 
 motion.mdl   Simulink library containing profile generators
              feedforward controllers and motion system models
             
Auxiliary and additional:
 make2.m      Matlab function to calculate timing for symmetrical 
              2nd order profiles plus profiles themselves
 make3.m      Matlab function to calculate timing for symmetrical 
              3d order profiles
 make4_it.m   Matlab function to calculate timing for symmetrical 
              4th order profiles using iterations to solve 
              3d order polynomial (obsolete)
 profile3.m   Matlab function to calculate symmetrical 3d order 
              profiles from timing info 

 makeplant.m  Matlab script to define several parameters and 
              state-space models for examples
              
 motion_example1.mdl   Simulink models giving examples of use
 motion_example2.mdl   of the motion library
 motion_example3.mdl
 ref4_rt.mdl           Real-time implementation for single axis
 ref4_xy.mdl           Real-time implementation for X-Y moves

***************************************************************
 
Installation and first use:
Extract files to work directory or any directory specified in 
matlabpath.
At Matlab command prompt: enter 'motion' to open motion library
and run examples.
For more details and command line use: see help make# and 
help profile#.
***************************************************************

Disclaimer:
This toolbox is provisional and experimental; hence no guarantees
whatsoever on results.
***************************************************************

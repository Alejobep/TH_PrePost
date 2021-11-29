# TH_PrePost
Scripts for pre- and post- processing numerical simulation files using TOUGH+HYDRATE v 1.5.

In addition to python, the only requirements are to have the pandas and numpy libraries installed.

Regarding the simulation files:

- The scripts currently handle only cartesian grids. 
- Input data file and the main output file need to have the same name, followed by the extensions "*.in" and "*.out", respectively.
- The input file needs to have the parameter MOP(19) set to 8. In this manner, it will produce an additional file containing the most important
properties printed in a format suitable for plotting.

As example, three different simulation cases are included to use for testing.

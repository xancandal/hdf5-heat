# Simulation of heat diffusion

Simulation which solves the heat diffusion in a volume, writing all intermediate results in a [HDF-5](http://www.hdfgroup.org/HDF5/) file, using a variable "temperature" with four dimensions (t, x, y, z). In addition, from the file HDF-5, the values of a sheet are extracted and written to a [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) file. 

Graphic representation of the time evolution of the temperature of a sheet extracted from the HDF-5 and NetCDF file is drawn. Note that it necessary to use a very small time differential to be stable (about 1e-9).

![sheet_2D](https://raw.githubusercontent.com/xancandal/hdf5-heat/master/Images/sheet_2D.png)

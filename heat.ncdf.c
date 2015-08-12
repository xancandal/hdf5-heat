#include <math.h>
#include <hdf5.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DX .01          // Interval size in x direction
#define DY .01          // Interval size in y direction
#define SECTION 50      // Value of z-plane slab
//// #define STEPS 500       // Timesteps to compute
#define STEPS 50    

// Filename hdf5
#define H5_FILE_NAME_HDF5 "data.h5"
#define ERR_HDF5 { H5Eprint( H5E_DEFAULT, stdout ); exit(2); }

// Filename NetCDF
#define NC_FILE_NAME_NETCDF "data.nc"
#define ERR_NETCDF(e) { printf( "Error: %s\n", nc_strerror(e) ); exit(2); }

int main() {

  double *data_out;       // buffer
  int nx = 1/DX;          // hyperslab and output buffer dimensions 
  int ny = 1/DY;
  // int i, j;
 
  // Allocate memory for data
  if((data_out = (double *)malloc(STEPS * nx * ny * sizeof(double))) == NULL)
    printf("Error malloc matrix data_out[%d]\n",nx * ny);

  // Create NetCDF file. NC_CLOBBER tells NetCDF to overwrite this file, if it already exists
  int ncid, retval;
  // if( retval = nc_create( NC_FILE_NAME_NETCDF, NC_CLOBBER|NC_NETCDF4, &ncid ) ) {
  if( retval = nc_create( NC_FILE_NAME_NETCDF, NC_CLOBBER, &ncid ) ) {
    ERR_NETCDF(retval);
  }

  // Define the x and y dimensions. NetCDF will hand back and ID for each.
  int x_dimid, y_dimid, t_dimid;
  if( retval = nc_def_dim( ncid, "x", nx, &x_dimid ) ) {
    ERR_NETCDF(retval);
  }
  if( retval = nc_def_dim( ncid, "y", ny, &y_dimid ) ) {
    ERR_NETCDF(retval);
  }
  
  // Define the t dimension at NetCDF.
  if( retval = nc_def_dim( ncid, "t", NC_UNLIMITED, &t_dimid ) ) {
    ERR_NETCDF(retval);
  }
  
  // Define coordinate variables for x and y at NetCDF
  int x_varid, y_varid, t_varid;
  if( retval = nc_def_var( ncid, "x", NC_DOUBLE, 1, &x_dimid, &x_varid ) ) {
    ERR_NETCDF(retval);
  }
  if( retval = nc_def_var( ncid, "y", NC_DOUBLE, 1, &y_dimid, &y_varid ) ) {
    ERR_NETCDF(retval);
  }
  if( retval = nc_def_var( ncid, "t", NC_DOUBLE, 1, &t_dimid, &t_varid ) ) {
    ERR_NETCDF(retval);
  }
  
  // Define the nc-variable to store temperature data. The t dimension should be the one which varies more slowly at NetCDF.
  int varid;
  int dimids[3] = { t_dimid, x_dimid, y_dimid };
  if( retval = nc_def_var( ncid, "temperature", NC_DOUBLE, 3, dimids, &varid )){
    ERR_NETCDF(retval);
  }

  // Write x, y, t and temperature units at NetCDF
  char * space_units = "meters";
  char * time_units = "seconds since start of the experiment";
  char * temp_units = "kelvin";
  if( retval = nc_put_att_text( ncid, x_varid, "units", strlen(space_units), space_units ) ) {
    ERR_NETCDF(retval);
  }
  if( retval = nc_put_att_text( ncid, y_varid, "units", strlen(space_units), space_units ) ) {
    ERR_NETCDF(retval);
  }
  if( retval = nc_put_att_text( ncid, t_varid, "units", strlen(time_units), time_units ) ) {
    ERR_NETCDF(retval);
  }
  if( retval = nc_put_att_text( ncid, varid, "units", strlen(temp_units), temp_units ) ) {
    ERR_NETCDF(retval);
  }
  double scale_factor = 300.0;
  if( retval = nc_put_att_double( ncid, varid, "scale_factor", NC_DOUBLE, 1, &scale_factor ) ) {
    ERR_NETCDF(retval);
  }

  // End define mode: this tells NetCDF that we are done defining metadata at NetCDF
  if( retval = nc_enddef( ncid ) ) {
    ERR_NETCDF(retval);
  }

  // Write x coordinates at NetCDF
  size_t pos;
  for( pos = 0; pos < nx; ++pos ) {
    double x = DX*pos;
    if( retval = nc_put_var1_double( ncid, x_varid, &pos, &x ) ) {
      ERR_NETCDF(retval);
    }
  }
  
  // Write y coordinates at NetCDF
  for( pos = 0; pos < ny; ++pos ) {
    double y = DY*pos;
    if( retval = nc_put_var1_double( ncid, y_varid, &pos, &y ) ) {
      ERR_NETCDF(retval);
    }
  }

  // Open an existing HDF5 file for Output buffer
  hid_t file_id = H5Fopen(H5_FILE_NAME_HDF5, H5F_ACC_RDONLY, H5P_DEFAULT);
  if( file_id < 0 ) { ERR_HDF5; }

  // Open an existing HDF5 dataset
  hid_t tempD = H5Dopen(file_id, "temperature", H5P_DEFAULT);
  if( tempD < 0 ) { ERR_HDF5; }

  // Returns an identifier for a copy of the dataspace for a dataset HDF5
  hid_t tempSel  = H5Dget_space (tempD);    /* dataspace handle */
  if( tempSel < 0 ) { ERR_HDF5; }

  // Returns the number of dimensions in the HDF5 dataspace if successful; otherwise returns a negative value
  int rank;
  rank  = H5Sget_simple_extent_ndims (tempSel);

  // Retrieves dataspace dimension size and maximum size HDF5
  hsize_t dims_out[2];    // HDF5 dataset dimensions
  hid_t status_n  = H5Sget_simple_extent_dims (tempSel, dims_out, NULL);
  if( status_n < 0 ) { ERR_HDF5; }
  
  // Display the number of dimensions in the HDF5 dataspace and the dataspace dimension size and maximum size 
  printf("\nRank: %d\nDimensions: %lu x %lu \n", rank, (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

  // Define hyperslab in the dataset HDF5
  hsize_t sel_offset_in[4] = {0,0,0,SECTION}; // The temperature value for the last interaction (STEPS) is chosen
  hsize_t sel_length_in[4] = {STEPS, nx, ny, 1};   
  H5Sselect_hyperslab( tempSel, H5S_SELECT_SET, sel_offset_in, NULL, sel_length_in, NULL );
  if( tempSel < 0 ) { ERR_HDF5; }

  // Define the memory dataspace HDF5
  hsize_t memSdim[4]={STEPS,nx,ny};
  hid_t memS = H5Screate_simple( 3, memSdim, NULL );
  if( memS < 0 ) { ERR_HDF5; }

  // Define memory HDF5 hyperslab
  hsize_t sel_offset_out[3] = {0,0,0};
  hsize_t sel_length_out[3] = {STEPS, nx, ny};
  H5Sselect_hyperslab( memS, H5S_SELECT_SET, sel_offset_out, NULL, sel_length_out, NULL );
  if( memS < 0 ) { ERR_HDF5; }

  // Read dataset tempD data from HDF5 hyperslab in the file into the hyperslab in memory
  hsize_t status = H5Dread (tempD, H5T_NATIVE_DOUBLE, memS, tempSel, H5P_DEFAULT, data_out);
  if( status < 0 ) { ERR_HDF5; }
  
  // printf ("Data:\n ");
  // for( i = 0; i < nx; ++i ) {
  //   for( j = 0; j < ny; ++j ) {
  //     printf("%f ", data_out[i*ny+j]);
  //   }
  //   printf("\n ");
  // }
  // printf("\n");

  // Write the data to the NETCDF file
  size_t corner_vector[3] = {0,0,0};
  size_t edge_lengths[3] = {STEPS, nx, ny};
  if(retval = nc_put_vara_double(ncid, varid, corner_vector, edge_lengths, data_out)){
    ERR_NETCDF(retval);
  }
  pos = 0;
  double tval = 0;
  if( retval = nc_put_var1_double( ncid, t_varid, &pos, &tval ) ) {
    ERR_NETCDF(retval);
  }

  // Close the HDF5 memspace
  if( H5Sclose( memS ) <  0 ) { ERR_HDF5; }

  // Close the HDF5 dataspace
  if( H5Sclose( tempSel ) <  0 ) { ERR_HDF5; }

  // Close the HDF5 dataset
  if( H5Dclose( tempD ) <  0 ) { ERR_HDF5; }

  // Close the HDF5 file
  if( H5Fclose( file_id ) <  0 ) { ERR_HDF5; }

  // Close the file. This frees up any internal NetCDF resources associated with the file, and flushes any buffers
  if( retval = nc_close( ncid ) ) {
    ERR_NETCDF(retval);
  }

  // Free memory
  free(data_out);

  return 0;
  
}



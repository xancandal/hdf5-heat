#include <math.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DX .01      // Interval size in x direction
#define DY .01      // Interval size in y direction
#define DZ .01      // Interval size in y direction
#define A 0.5       // Diffusion constant
//// #define STEPS 500   // Timesteps to compute
#define STEPS 50  

#define H5_FILE_NAME "data.h5"
#define ERR { H5Eprint( H5E_DEFAULT, stdout ); exit(2); }

// Allocate memory function
void reserva_memoria( double **u, int nx, int ny, int nz ) {
    if((*u = (double *)malloc(nx * ny * nz * sizeof(double))) == NULL)
    printf("Error malloc matrix u[%d]\n",nx * ny * nz);
}

int main() {
  double *u, *ui;
  int i,j,k,t;

  int nx = 1/DX;  // Domain
  int ny = 1/DY;
  int nz = 1/DZ;

  // To save CPU cycles, we'll compute Delta x^2, Delta y^2 and Delta z^2 only once
  double dx2 = DX*DX;
  double dy2 = DY*DY;
  double dz2 = DZ*DZ;
  //double dt = dx2*dy2*dz2/(2*A*(dx2+dy2+dz2));
  double dt = 0.00003;

 // Allocate memory
  reserva_memoria(&u,nx,ny,nz);
  reserva_memoria(&ui,nx,ny,nz);

  // Initialize matrices, Superficie o Volumen
  for( i = 0; i < nx; ++i ) {
    for( j = 0; j < ny; ++j ) {
      for( k = 0; k < nz; ++k ) {
        double dist_to_center = pow(i*DX-0.5,2)+pow(j*DY-0.5,2)+pow(k*DZ-0.5,2);
        if( ( dist_to_center < 0.1 ) && ( dist_to_center > 0.05 ) ) {
          ui[i*ny*nz+j*nz+k] = 1.0;
        }
        else {
          ui[i*ny*nz+j*nz+k] = 0.0;
        }
      }
    }
  }
  
  // Create HDF5 file. The H5F_ACC_TRUNC parameter tells HDF5 to overwrite the this file if it already exists.
  hid_t file_id;
  file_id = H5Fcreate( H5_FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  if( file_id < 0 ) { ERR; }
  
  // Create a dataspace for the temperature. This time around we need an unlimited time dimension.
  hsize_t curr_dims[4] = { 0, nx, ny, nz };
  hsize_t max_dims[4] = {H5S_UNLIMITED, nx, ny, nz};
  hid_t tempS = H5Screate_simple( 4, curr_dims, max_dims );
  if( tempS < 0 ) { ERR; }

  // Create the Dimension Scale for the X variable and write values
  double xds_data[nx];
  hsize_t xdims = nx;
  hid_t xS = H5Screate_simple( 1, &xdims, NULL );
  hid_t xD = H5Dcreate( file_id, "x", H5T_IEEE_F64BE, xS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  for( i = 0; i < nx; ++i ) {
    xds_data[i] = DX*i;
  }
  herr_t status = H5Dwrite( xD, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xds_data);
  status = H5DSset_scale( xD, "x" );

  // Create the Dimension Scale for the Y variable and write values
  double yds_data[ny];
  hsize_t ydims = ny;
  hid_t yS = H5Screate_simple( 1, &ydims, NULL );
  hid_t yD = H5Dcreate( file_id, "y", H5T_IEEE_F64BE, yS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  for( i = 0; i < ny; ++i ) {
    yds_data[i] = DY*i;
  }
  status = H5Dwrite( yD, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yds_data);
  status = H5DSset_scale( yD, "y" );

  // Create the Dimension Scale for the Z variable and write values
  double zds_data[nz];
  hsize_t zdims = nz;
  hid_t zS = H5Screate_simple( 1, &zdims, NULL );
  hid_t zD = H5Dcreate( file_id, "z", H5T_IEEE_F64BE, zS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  for( i = 0; i < nz; ++i ) {
    zds_data[i] = DZ*i;
  }
  status = H5Dwrite( zD, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zds_data);
  status = H5DSset_scale( zD, "z" );

  // Create the Dimension Scale for time
  hid_t tS = H5Screate_simple( 1, curr_dims, max_dims );
  hid_t tCPL = H5Pcreate( H5P_DATASET_CREATE );
  hsize_t tchunk_dims = 1;
  status = H5Pset_chunk( tCPL, 1, &tchunk_dims );
  hid_t tD = H5Dcreate( file_id, "t", H5T_IEEE_F64BE, tS, H5P_DEFAULT, tCPL,H5P_DEFAULT );
  status = H5DSset_scale( tD, "t" );

  // Create the temperature dataset into the root group of the file. Given that
  // the dataspace for this dataset has unlimited dimensions, it is mandatory
  // that the dataset is chunked.
  hid_t tempCPL = H5Pcreate( H5P_DATASET_CREATE );
  hsize_t chunk_dims[4] = { 1, nx, ny, nz };
  H5Pset_chunk( tempCPL, 4, chunk_dims );

  // Compress data temperature
  // H5Pset_deflate(tempCPL, 9);

  hid_t tempD = H5Dcreate( file_id, "temperature", H5T_IEEE_F64BE, tempS,H5P_DEFAULT, tempCPL, H5P_DEFAULT );
  if( tempD < 0 ) { ERR; }
  H5DSattach_scale( tempD, tD, 0 );
  H5DSattach_scale( tempD, xD, 1 );
  H5DSattach_scale( tempD, yD, 2 );
  H5DSattach_scale( tempD, zD, 3 );
  
  // Write the initial state
  ++curr_dims[0];
  H5Dset_extent( tempD, curr_dims );

  hid_t tempSel = H5Dget_space( tempD );
  if( tempSel < 0 ) { ERR; }

  hsize_t sel_offset[4] = {0,0,0,0};
  hsize_t sel_length[4] = {1, nx, ny, nz};
  hsize_t memSdim[3]={nx,ny,nz};
  hid_t memS = H5Screate_simple( 3, memSdim, NULL );
  if( memS < 0 ) { ERR; }
  H5Sselect_hyperslab( tempSel, H5S_SELECT_SET, sel_offset, NULL, sel_length, NULL );
  if( tempSel < 0 ) { ERR; }
  H5Dwrite( tempD, H5T_NATIVE_DOUBLE, memS, tempSel, H5P_DEFAULT, ui );

  H5Dset_extent( tD, curr_dims );
  
  // Write also to the time dimension scale
  hid_t tSel = H5Dget_space( tD );
  H5Sselect_hyperslab(tSel, H5S_SELECT_SET, sel_offset, NULL, sel_length, NULL);
  double tval = 0;
  hid_t scalarS = H5Screate( H5S_SCALAR );
  H5Dwrite( tD, H5T_NATIVE_DOUBLE, scalarS, tSel, H5P_DEFAULT, &tval );

  // Calculate all timesteps
  //// double du = 1;
  //// t = 0;
  //// while( du > 0.1 ) {
  for( t = 0; t < STEPS; ++t ) {
    printf( "TIMESTEP = %d\n", t );
    //// du = 0;
    for( i = 1; i < nx-1; ++i ) {
      for( j = 1; j < ny-1; ++j ) {
        for( k = 1; k < nz-1; ++k ) {
          u[i*ny*nz+j*nz+k] = ui[i*ny*nz+j*nz+k] + A*dt*( 
             ( ui[(i+1)*ny*nz+j*nz+k] - 2*ui[i*ny*nz+j*nz+k] + ui[(i-1)*ny*nz+j*nz+k] )/dx2
           + ( ui[i*ny*nz+(j+1)*nz+k] - 2*ui[i*ny*nz+j*nz+k] + ui[i*ny*nz+(j-1)*nz+k] )/dy2 
           + ( ui[i*ny*nz+j*nz+(k+1)] - 2*ui[i*ny*nz+j*nz+k] + ui[i*ny*nz+j*nz+(k-1)] )/dz2 );
        //// du += fabs( u[i*ny*nz+j*nz+k] - ui[i*ny*nz+j*nz+k] );
        }
      }
    }
    
    //// printf( "TIMESTEP = %d (du=%lf)\n", t, du );
    printf( "TIMESTEP = %d \n", t );
    memcpy( ui, u, nx*ny*nz*sizeof(double) );
    
    // Write this t point to the dimension scale
    ++curr_dims[0];
    H5Dset_extent( tD, curr_dims );
    tSel = H5Dget_space( tD );
    if( tSel < 0 ) { ERR; }
    hsize_t offset = t+1;
    status = H5Sselect_hyperslab( tSel, H5S_SELECT_SET, &offset, NULL, sel_length, NULL );
    if( tSel < 0 ) { ERR; }
    tval = (t+1)*dt;
    status = H5Dwrite( tD, H5T_NATIVE_DOUBLE, scalarS, tSel, H5P_DEFAULT, &tval );
    if( status < 0 ) { ERR; }
    
    // Write the temperature data to the file
    H5Dset_extent( tempD, curr_dims );
    tempSel = H5Dget_space( tempD ); // The dataspace has changed!
    if( tempSel < 0 ) { ERR; }
    sel_offset[0] = t+1;
    H5Sselect_hyperslab( tempSel, H5S_SELECT_SET, sel_offset, NULL, sel_length, NULL );
    if( tempSel < 0 ) { ERR; }
    status = H5Dwrite( tempD, H5T_NATIVE_DOUBLE, memS, tempSel, H5P_DEFAULT, u);
    if( status < 0 ) { ERR; }
    
    //// ++t;

  }
  
  // Close the file. This frees up any internal NetCDF resources associated
  // with the file, and flushes any buffers.
  if( H5Fclose( file_id ) <  0 ) { ERR; }
     
  // Free memory
  free(u);
  free(ui);

  return 0;
}

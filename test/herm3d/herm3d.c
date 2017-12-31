/* Wind simulation library.
 * Copyright (C) 2013  Michael Andre.
 * Copyright (C) 2013  Technische Universitaet Muenchen.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "herm.h"
#include "rv.h"

#include <stdio.h>
#include <math.h>

void copy (complex_t * in, complex_t * out, dft_ptr_t size)
{
  dft_ptr_t i;
  
  for (i=0; i < size; i++) {
    out[i][0] = in[i][0];
    out[i][1] = in[i][1];
  }
}

double normrc (complex_t * fc1, complex_t * fc2, dft_ptr_t size)
{
  dft_ptr_t i;
  double y, ymax;

  ymax = 0.;
  for (i=0; i < size; i++) {
    y = (fc1[i][0] - fc2[i][0])*(fc1[i][0] - fc2[i][0]) + (fc1[i][1] - fc2[i][1])*(fc1[i][1] - fc2[i][1]);
    ymax = (ymax > y) ? ymax : y;
  }
  return sqrt(ymax);
}

void normalize (complex_t * a, dft_ptr_t size, dft_ptr_t n)
{
  double rfac = 1. / sqrt((double) n);
  dft_ptr_t i;

  for (i=0; i < size; i++) {
    a[i][0] = rfac * a[i][0];
    a[i][1] = rfac * a[i][1];
  }
}

void apply_symmetry(complex_t * fc, dft_ptr_t nx, dft_ptr_t ny, dft_ptr_t nz, dft_ptr_t local_x_start, dft_ptr_t local_nx)
{
  dft_ptr_t i, j, k, l, iptr, jptr, ptr;
  complex_t * slice;
  herm_t * herm;

  herm = herm_init(nx, ny, local_x_start, local_nx);
  slice = herm->slice;

  // set Hermitian planes k = 0 
  k = 0;
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        slice[i*ny + j][l] = fc[ptr][l];
      }
    }
  }

  // apply Hermitian symmetry k = 0 
  herm_execute(herm);

  // copy Hermitian planes k = 0                          
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        fc[ptr][l] = slice[i*ny + j][l];
      }
    }
  }

  // set Hermitian planes k = nz/2      
  k = nz/2;
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        slice[i*ny + j][l] = fc[ptr][l];
      }
    }  
  }  

  // apply Hermitian symmetry k = nz/2
  herm_execute(herm);

  // copy Hermitian planes k = nz/2     
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        fc[ptr][l] = slice[i*ny + j][l];
      }
    }  
  }  

  herm_destroy(herm);
}

int main(int argc, char **argv)
{
  if (argc != 4) {
    printf("Wrong number of arguments\n");
    return 1;
  }

  int rank;
  double lnorm, gnorm;
  dft_ptr_t i, j, k, size;
  dft_ptr_t rx, ry, rz, nx, ny, nz, local_x_start, local_nx;
  void * cin, * rcinout;
  fftw_plan planr2c, planc2r;

  rx = atoi(argv[1]);
  ry = atoi(argv[2]);
  rz = atoi(argv[3]);

  nx = 1;
  for (i=0; i < rx; i++)
    nx = 2 * nx;

  ny = 1;
  for (i=0; i < ry; i++)
    ny = 2 * ny;

  nz = 1;
  for (i=0; i < rz; i++)
    nz = 2 * nz;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fftw_mpi_init();
  size = fftw_mpi_local_size_3d(nx, ny, nz/2 + 1, MPI_COMM_WORLD, &local_nx, &local_x_start);
#else /* DON'T HAVE_MPI */
  rank = 0;
  local_nx = nx;
  local_x_start = 0;
  size = local_nx * ny * (nz/2 + 1);
#endif /* HAVE_MPI */

  rv_init(rank);

  cin  = fftw_malloc( sizeof(complex_t) * size );
  rcinout = fftw_malloc( sizeof(complex_t) * size );
  size = local_nx * ny * (nz/2 + 1);

  // generate complex numbers
  if (rank == 0) printf("Generating complex random numbers\n");
  rv_normal(cin, 2 * size);

  if (rank == 0) printf("Applying Hermitian symmetry\n");
  apply_symmetry(cin, nx, ny, nz, local_x_start, local_nx);

  copy(cin, rcinout, size);

  // inverse dft
#ifdef HAVE_MPI
  planc2r = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, rcinout, rcinout, MPI_COMM_WORLD, FFTW_ESTIMATE);
#else /* DON'T HAVE_MPI */
  planc2r = fftw_plan_dft_c2r_3d(nx, ny, nz, rcinout, rcinout, FFTW_ESTIMATE);
#endif /* HAVE_MPI */

  if (rank == 0) printf("Calculating inverse DFT\n");
  fftw_execute(planc2r);

  // forward dft
#ifdef HAVE_MPI
  planr2c = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, rcinout, rcinout, MPI_COMM_WORLD, FFTW_ESTIMATE);
#else /* DON'T HAVE_MPI */
  planr2c = fftw_plan_dft_r2c_3d(nx, ny, nz, rcinout, rcinout, FFTW_ESTIMATE);
#endif /* HAVE_MPI */

  if (rank == 0) printf("Calculating forward DFT\n");
  fftw_execute(planr2c);
  normalize(rcinout, size, nx*ny*nz); // backward
  normalize(rcinout, size, nx*ny*nz); // forward
  
  lnorm = normrc(cin, rcinout, size);
#ifdef HAVE_MPI
  MPI_Reduce(&lnorm, &gnorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else /* DON'T HAVE_MPI */
  gnorm = lnorm;
#endif /* HAVE_MPI */

  if (rank == 0) printf("max||F^-1F(f) - f|| = %e\n", gnorm);

  fftw_destroy_plan(planr2c);
  fftw_destroy_plan(planc2r);

#ifdef HAVE_MPI
  fftw_mpi_cleanup();
  MPI_Finalize();
#endif /* HAVE_MPI */

  fftw_free(cin);
  fftw_free(rcinout);

  return 0;
}

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

double normc (complex_t * a, dft_ptr_t size)
{
  dft_ptr_t i;
  double y, ymax;

  ymax = 0.;
  for (i=0; i < size; i++) {
    y = a[i][1] * a[i][1];
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

void copy (double * in, double * out, dft_ptr_t size)
{
  dft_ptr_t i;
  for (i=0; i < size; i++) {
    out[i] = in[i];
  }
}

int main(int argc, char **argv)
{
  if (argc != 3) {
    printf("Wrong number of arguments\n");
    return 1;
  }

  herm_t * herm;
  int rank;
  double lnorm, gnorm;
  dft_ptr_t i, j, size;
  dft_ptr_t rx, ry, nx, ny, local_x_start, local_nx;
  void * cinout;
  fftw_plan plan;

  rx = atoi(argv[1]);
  ry = atoi(argv[2]);

  nx = 1;
  for (i=0; i < rx; i++)
    nx = 2 * nx;

  ny = 1;
  for (i=0; i < ry; i++)
    ny = 2 * ny;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fftw_mpi_init();
  size = fftw_mpi_local_size_2d(nx, ny, MPI_COMM_WORLD, &local_nx, &local_x_start);
#else /* DON'T HAVE_MPI */
  rank = 0;
  local_nx = nx;
  local_x_start = 0;
  size = local_nx * ny;
#endif /* HAVE_MPI */

  rv_init(rank);

  cinout = fftw_malloc( sizeof(complex_t) * size );
  size = local_nx * ny;

  // generate complex numbers
  rv_normal(cinout, 2 * size);

  lnorm = normc(cinout, size);
#ifdef HAVE_MPI
  MPI_Reduce(&lnorm, &gnorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else /* DON'T HAVE_MPI */
  gnorm = lnorm;
#endif /* HAVE_MPI */

  if (rank == 0) printf("Imaginary norm of raw data = %e\n", gnorm);
  if (rank == 0) printf("Applying Hermitian symmetry\n");

  // apply Hermitian symmetry
  herm = herm_init(nx, ny, local_x_start, local_nx);

  copy(cinout, herm->slice, 2 * size);
  herm_execute(herm);
  copy(herm->slice, cinout, 2 * size);
  herm_destroy(herm);

  // inverse dft
#ifdef HAVE_MPI
  plan = fftw_mpi_plan_dft_2d(nx, ny, cinout, cinout, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
#else /* DON'T HAVE_MPI */
  plan = fftw_plan_dft_2d(nx, ny, cinout, cinout, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif /* HAVE_MPI */

  if (rank == 0) printf("Calculating inverse DFT\n");

  fftw_execute(plan);
  normalize(cinout, size, nx*ny);

  lnorm = normc(cinout, size);
#ifdef HAVE_MPI
  MPI_Reduce(&lnorm, &gnorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#else /* DON'T HAVE_MPI */
  gnorm = lnorm;
#endif /* HAVE_MPI */

  if (rank == 0) printf("Imaginary norm of inverse DFT = %e\n", gnorm);

  fftw_destroy_plan(plan);

#ifdef HAVE_MPI
  fftw_mpi_cleanup();
  MPI_Finalize();
#endif /* HAVE_MPI */

  fftw_free(cinout);

  return 0;
}

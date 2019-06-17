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

#include "dft.h"

void dft_init_mpi(void)
{
#ifdef HAVE_MPI
  fftw_mpi_init();
#endif /* HAVE_MPI */
}

void dft_finalize_mpi(void)
{
#ifdef HAVE_MPI
  fftw_mpi_cleanup();
#endif /* HAVE_MPI */
}

/**
 * dft_init3D:
 * Initialize a 3-dimensional field for a complex-to-real inverse DFT.
 */
dft_t * dft_init3D(dft_ptr_t nx, dft_ptr_t ny, dft_ptr_t nz)
{
  dft_ptr_t size;
  dft_t * dft;

  dft = (dft_t *) malloc( sizeof(dft_t) );

#ifdef HAVE_MPI
  size = fftw_mpi_local_size_3d(nx, ny, nz/2 + 1, MPI_COMM_WORLD, &(dft->local_nx), &(dft->local_x_start));
#else /* DON'T HAVE_MPI */
  dft->local_nx = nx;
  dft->local_x_start = 0;
  size = nx * ny * (nz/2 + 1);
#endif /* HAVE_MPI */

  dft->field = fftw_malloc( 2 * sizeof(double) * size );

#ifdef HAVE_MPI
  dft->plan = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, dft->field, dft->field, MPI_COMM_WORLD, FFTW_ESTIMATE);
#else /* DON'T HAVE_MPI */
  dft->plan = fftw_plan_dft_c2r_3d(nx, ny, nz, dft->field, dft->field, FFTW_ESTIMATE);
#endif /* HAVE_MPI */

  return dft;
}

/**
 * dft_init2D:
 * Initialize a 2-dimensional field for a complex-to-real inverse DFT.
 */
dft_t * dft_init2D(dft_ptr_t nx, dft_ptr_t nz)
{
  dft_ptr_t size;
  dft_t * dft;

  dft = (dft_t *) malloc( sizeof(dft_t) );

#ifdef HAVE_MPI
  size = fftw_mpi_local_size_2d(nx, nz/2 + 1, MPI_COMM_WORLD, &(dft->local_nx), &(dft->local_x_start));
#else /* DON'T HAVE_MPI */
  dft->local_nx = nx;
  dft->local_x_start = 0;
  size = nx * (nz/2 + 1);
#endif /* HAVE_MPI */

  dft->field = fftw_malloc( 2 * sizeof(double) * size );

#ifdef HAVE_MPI
  dft->plan = fftw_mpi_plan_dft_c2r_2d(nx, nz, dft->field, dft->field, MPI_COMM_WORLD, FFTW_ESTIMATE);
#else /* DON'T HAVE_MPI */
  dft->plan = fftw_plan_dft_c2r_2d(nx, nz, dft->field, dft->field, FFTW_ESTIMATE);
#endif /* HAVE_MPI */

  return dft;
}

void dft_execute(dft_t * dft)
{
  fftw_execute(dft->plan);
}

void dft_destroy(dft_t * dft)
{
  fftw_destroy_plan(dft->plan);
  fftw_free(dft->field);
  free(dft);
}

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

/**
 * herm_init:
 * Initialize a 2-dimensional, periodic field on which Hermitian symmetry will be enforced.
 */
herm_t * herm_init(dft_ptr_t nx, dft_ptr_t ny, dft_ptr_t local_x_start, dft_ptr_t local_nx)
{
  herm_t * herm;
  int p;
  dft_ptr_t slice_out_nx, slice_in_nx, size;

  herm = (herm_t *) malloc( sizeof(herm_t) );

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &(herm->rank));
  MPI_Comm_size(MPI_COMM_WORLD, &(herm->nproc));
#else /* DON'T HAVE_MPI */
  herm->rank = 0;
  herm->nproc = 1;
#endif /* HAVE_MPI */

  herm->slice = NULL;
  herm->slice_in = NULL;
  herm->slice_out = NULL;
  herm->nx = nx;
  herm->ny = ny;
  herm->local_x_start     = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_nx          = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_out_floor = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_out_ceil  = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_in_floor  = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_in_ceil   = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );

#ifdef HAVE_MPI
  if (sizeof(dft_ptr_t) == sizeof(int)) {
    MPI_Allgather(&local_x_start, 1, MPI_INT, herm->local_x_start, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_INT, herm->local_nx, 1, MPI_INT, MPI_COMM_WORLD);
  }
  else if (sizeof(dft_ptr_t) == sizeof(long)) {
    MPI_Allgather(&local_x_start, 1, MPI_LONG, herm->local_x_start, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_LONG, herm->local_nx, 1, MPI_LONG, MPI_COMM_WORLD);
  }
  else if (sizeof(dft_ptr_t) == sizeof(long long)) {
    MPI_Allgather(&local_x_start, 1, MPI_LONG_LONG, herm->local_x_start, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_LONG_LONG, herm->local_nx, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
  }
  else {
    herm_destroy(herm);
    return NULL;
  }
#else /* DON'T HAVE_MPI */
  herm->local_x_start[0] = local_x_start;
  herm->local_nx[0] = local_nx;
#endif /* HAVE_MPI */

  for (p=0; p < herm->nproc; p++) {
    herm->local_x_out_floor[p] = herm->local_x_start[p];
    herm->local_x_out_ceil[p] = min(nx/2, herm->local_x_start[p]+herm->local_nx[p]);
    herm->local_x_in_floor[p] = max(nx/2+1, herm->local_x_start[p]);
    herm->local_x_in_ceil[p] = herm->local_x_start[p]+herm->local_nx[p];
  }
  herm->local_x_out_floor[0] = 1;

  slice_out_nx = herm->local_x_out_ceil[herm->rank] - herm->local_x_out_floor[herm->rank];
  slice_in_nx = herm->local_x_in_ceil[herm->rank] - herm->local_x_in_floor[herm->rank];

  size = herm->local_nx[herm->rank] * herm->ny;
  if (size > 0)
    herm->slice = fftw_malloc( sizeof(complex_t) * size );

  size = slice_out_nx * herm->ny;
  if (size > 0)
    herm->slice_out = fftw_malloc( sizeof(complex_t) * size );

  size = slice_in_nx * herm->ny;
  if (size > 0)
    herm->slice_in = fftw_malloc( sizeof(complex_t) * size );

  return herm;
}

/**
 * herm_init2D:
 * Initialize a one-dimensional, periodic field on which Hermitian symmetry will be enforced.
 */
herm_t * herm_init2D(dft_ptr_t nx, dft_ptr_t local_x_start, dft_ptr_t local_nx)
{
  herm_t * herm;
  int p;
  dft_ptr_t slice_out_nx, slice_in_nx, size;

  herm = (herm_t *) malloc( sizeof(herm_t) );

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &(herm->rank));
  MPI_Comm_size(MPI_COMM_WORLD, &(herm->nproc));
#else /* DON'T HAVE_MPI */
  herm->rank = 0;
  herm->nproc = 1;
#endif /* HAVE_MPI */

  herm->slice = NULL;
  herm->slice_in = NULL;
  herm->slice_out = NULL;
  herm->nx = nx;
  herm->local_x_start     = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_nx          = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_out_floor = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_out_ceil  = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_in_floor  = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );
  herm->local_x_in_ceil   = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * herm->nproc );

#ifdef HAVE_MPI
  if (sizeof(dft_ptr_t) == sizeof(int)) {
    MPI_Allgather(&local_x_start, 1, MPI_INT, herm->local_x_start, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_INT, herm->local_nx, 1, MPI_INT, MPI_COMM_WORLD);
  }
  else if (sizeof(dft_ptr_t) == sizeof(long)) {
    MPI_Allgather(&local_x_start, 1, MPI_LONG, herm->local_x_start, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_LONG, herm->local_nx, 1, MPI_LONG, MPI_COMM_WORLD);
  }
  else if (sizeof(dft_ptr_t) == sizeof(long long)) {
    MPI_Allgather(&local_x_start, 1, MPI_LONG_LONG, herm->local_x_start, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_LONG_LONG, herm->local_nx, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
  }
  else {
    herm_destroy(herm);
    return NULL;
  }
#else /* DON'T HAVE_MPI */
  herm->local_x_start[0] = local_x_start;
  herm->local_nx[0] = local_nx;
#endif /* HAVE_MPI */

  for (p=0; p < herm->nproc; p++) {
    herm->local_x_out_floor[p] = herm->local_x_start[p];
    herm->local_x_out_ceil[p] = min(nx/2, herm->local_x_start[p]+herm->local_nx[p]);
    herm->local_x_in_floor[p] = max(nx/2+1, herm->local_x_start[p]);
    herm->local_x_in_ceil[p] = herm->local_x_start[p]+herm->local_nx[p];
  }
  herm->local_x_out_floor[0] = 1;

  slice_out_nx = herm->local_x_out_ceil[herm->rank] - herm->local_x_out_floor[herm->rank];
  slice_in_nx = herm->local_x_in_ceil[herm->rank] - herm->local_x_in_floor[herm->rank];

  size = herm->local_nx[herm->rank];
  if (size > 0)
    herm->slice = fftw_malloc( sizeof(complex_t) * size );

  size = slice_out_nx;
  if (size > 0)
    herm->slice_out = fftw_malloc( sizeof(complex_t) * size );

  size = slice_in_nx;
  if (size > 0)
    herm->slice_in = fftw_malloc( sizeof(complex_t) * size );

  return herm;
}

/**
 * herm_execute:
 * Apply Hermitian symmetry.
 */
void herm_execute(herm_t * herm)
{
  double *   slice               = herm->slice;
  double *   slice_in            = herm->slice_in;
  double *   slice_out           = herm->slice_out;
  dft_ptr_t  nx                  = herm->nx;
  dft_ptr_t  ny                  = herm->ny;
  dft_ptr_t  local_x_start       = herm->local_x_start[herm->rank];
  dft_ptr_t  local_x_end         = herm->local_x_start[herm->rank] + herm->local_nx[herm->rank];
  dft_ptr_t  slice_out_nx        = herm->local_x_out_ceil[herm->rank] - herm->local_x_out_floor[herm->rank];
  dft_ptr_t  slice_in_nx         = herm->local_x_in_ceil[herm->rank] - herm->local_x_in_floor[herm->rank];
  dft_ptr_t  slice_out_x_start   = herm->local_x_out_floor[herm->rank];
  dft_ptr_t  slice_in_x_start    = herm->local_x_in_floor[herm->rank];
  int ps, pr;
  double * tmp;
  dft_ptr_t i, j, ptr1, ptr2, iptr, jptr, size;
  dft_ptr_t nx_overlap, local_x_out_start, local_x_in_start;


  // i = 0
  if (herm->rank == 0) { // perodic-x and Hermitian symmetry
    for (j=1; j < ny/2; j++) {
      jptr = 2*j;
      slice[2*ny - jptr] = slice[jptr]; // symmetric real
      slice[2*ny - jptr + 1] = -slice[jptr + 1]; // symmetric imag
    }
    slice[1] = 0.;
    slice[ny + 1] = 0.;
  }

  // i = nx/2
  if (nx/2 >= local_x_start && nx/2 < local_x_end) {
    iptr = 2 * (nx/2 - local_x_start) * ny;
    for (j=1; j < ny/2; j++) {
      jptr = 2*j;
      slice[iptr + 2*ny - jptr] = slice[iptr + jptr]; 
      slice[iptr + 2*ny - jptr + 1] = -slice[iptr + jptr + 1]; 
    }
    slice[iptr + 1] = 0.;
    slice[iptr + ny + 1] = 0.;
  }

  // fill slice_out
  ptr1 = 2 * (slice_out_x_start - local_x_start) * ny;
  size = 2 * slice_out_nx * ny;
  for (i=0; i < size; i++) {
    slice_out[i] = slice[ptr1 + i];
  }

  // apply local symmetry
  tmp = (double *) malloc( sizeof(complex_t) * slice_out_nx * ny );
  ptr1 = 2 * ( (slice_out_nx-1) * ny + ny );

  for (i=0; i < slice_out_nx; i++) {
    iptr = 2 * i * ny;
    tmp[ptr1 - 2*ny - iptr] = slice_out[iptr];
    tmp[ptr1 - 2*ny - iptr + 1] = -slice_out[iptr + 1];
    for (j=1; j < ny; j++) {
      jptr = 2 * j;
      tmp[ptr1 - iptr - jptr] = slice_out[iptr + jptr];
      tmp[ptr1 - iptr - jptr + 1] = -slice_out[iptr + jptr + 1];
    }
  }

  for (i=0; i < slice_out_nx; i++) {
    iptr = 2 * i * ny;
    for (j=0; j < ny; j++) {
      jptr = 2 * j;
      slice_out[iptr + jptr] = tmp[iptr + jptr];
      slice_out[iptr + jptr + 1] = tmp[iptr + jptr + 1];
    }
  }

  free(tmp);

  // apply global symmetry
  for (ps=0; ps < herm->nproc; ps++) {
    for (pr=0; pr < herm->nproc; pr++) {
      local_x_out_start =  herm->local_x_out_floor[ps] + max(herm->local_x_out_ceil[ps] - (nx+1 - herm->local_x_in_floor[pr]), 0);
      local_x_in_start = nx+1 - min(nx+1 - herm->local_x_in_floor[pr], herm->local_x_out_ceil[ps]);
      nx_overlap = nx+1 - local_x_in_start - max(nx+1 - herm->local_x_in_ceil[pr], herm->local_x_out_floor[ps]);
      ptr1 = 2 * (local_x_in_start - herm->local_x_in_floor[pr]) * ny;
      ptr2 = 2 * (local_x_out_start - herm->local_x_out_floor[ps]) * ny;
      size = 2 * nx_overlap * ny;
      
      if (nx_overlap > 0) {
#ifdef HAVE_MPI
	if (herm->rank == ps && herm->rank != pr) {
	  MPI_Send(slice_out+ptr2, size, MPI_DOUBLE, pr, ps, MPI_COMM_WORLD);
	}
	else if (herm->rank != ps && herm->rank == pr) {
	  MPI_Recv(slice_in+ptr1, size, MPI_DOUBLE, ps, ps, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if (herm->rank == ps && herm->rank == pr) {
	  for (i=0; i < size; i++) {
	    slice_in[ptr1 + i] = slice_out[ptr2 + i];
	  }
	}
#else /* DON'T HAVE_MPI */
	for (i=0; i < size; i++) {
	  slice_in[ptr1 + i] = slice_out[ptr2 + i];
	}
#endif /* HAVE_MPI */
      } /* end if */
    }
  }

  // copy slice_in
  ptr1 = 2 * (slice_in_x_start - local_x_start) * ny;
  size = 2 * slice_in_nx * ny;
  for (i=0; i < size; i++) {
    slice[ptr1 + i] = slice_in[i];
  }
}

/**
 * herm_execute2D:
 * Apply Hermitian symmetry.
 */
void herm_execute2D(herm_t * herm)
{
  double *   slice               = herm->slice;
  double *   slice_in            = herm->slice_in;
  double *   slice_out           = herm->slice_out;
  dft_ptr_t  nx                  = herm->nx;
  dft_ptr_t  local_x_start       = herm->local_x_start[herm->rank];
  dft_ptr_t  local_x_end         = herm->local_x_start[herm->rank] + herm->local_nx[herm->rank];
  dft_ptr_t  slice_out_nx        = herm->local_x_out_ceil[herm->rank] - herm->local_x_out_floor[herm->rank];
  dft_ptr_t  slice_in_nx         = herm->local_x_in_ceil[herm->rank] - herm->local_x_in_floor[herm->rank];
  dft_ptr_t  slice_out_x_start   = herm->local_x_out_floor[herm->rank];
  dft_ptr_t  slice_in_x_start    = herm->local_x_in_floor[herm->rank];
  int ps, pr;
  double * tmp;
  dft_ptr_t i, j, ptr1, ptr2, iptr, jptr, size;
  dft_ptr_t nx_overlap, local_x_out_start, local_x_in_start;


  // i = 0
  slice[1] = 0.;

  // i = nx/2
  if (nx/2 >= local_x_start && nx/2 < local_x_end) {
    iptr = 2 * (nx/2 - local_x_start);
    slice[iptr + 1] = 0.;
  }

  // fill slice_out
  ptr1 = 2 * (slice_out_x_start - local_x_start);
  size = 2 * slice_out_nx;
  for (i=0; i < size; i++) {
    slice_out[i] = slice[ptr1 + i];
  }

  // apply local symmetry
  tmp = (double *) malloc( sizeof(complex_t) * slice_out_nx );
  ptr1 = 2 * (slice_out_nx-1);

  for (i=0; i < slice_out_nx; i++) {
    iptr = 2 * i;
    tmp[ptr1 - iptr] = slice_out[iptr];
    tmp[ptr1 - iptr + 1] = -slice_out[iptr + 1];
  }

  for (i=0; i < slice_out_nx; i++) {
    iptr = 2 * i;
    slice_out[iptr] = tmp[iptr];
    slice_out[iptr + 1] = tmp[iptr + 1];
  }

  free(tmp);

  // apply global symmetry
  for (ps=0; ps < herm->nproc; ps++) {
    for (pr=0; pr < herm->nproc; pr++) {
      local_x_out_start =  herm->local_x_out_floor[ps] + max(herm->local_x_out_ceil[ps] - (nx+1 - herm->local_x_in_floor[pr]), 0);
      local_x_in_start = nx+1 - min(nx+1 - herm->local_x_in_floor[pr], herm->local_x_out_ceil[ps]);
      nx_overlap = nx+1 - local_x_in_start - max(nx+1 - herm->local_x_in_ceil[pr], herm->local_x_out_floor[ps]);
      ptr1 = 2 * (local_x_in_start - herm->local_x_in_floor[pr]);
      ptr2 = 2 * (local_x_out_start - herm->local_x_out_floor[ps]);
      size = 2 * nx_overlap;
      
      if (nx_overlap > 0) {
#ifdef HAVE_MPI
  if (herm->rank == ps && herm->rank != pr) {
    MPI_Send(slice_out+ptr2, size, MPI_DOUBLE, pr, ps, MPI_COMM_WORLD);
  }
  else if (herm->rank != ps && herm->rank == pr) {
    MPI_Recv(slice_in+ptr1, size, MPI_DOUBLE, ps, ps, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (herm->rank == ps && herm->rank == pr) {
    for (i=0; i < size; i++) {
      slice_in[ptr1 + i] = slice_out[ptr2 + i];
    }
  }
#else /* DON'T HAVE_MPI */
  for (i=0; i < size; i++) {
    slice_in[ptr1 + i] = slice_out[ptr2 + i];
  }
#endif /* HAVE_MPI */
      } /* end if */
    }
  }

  // copy slice_in
  ptr1 = 2 * (slice_in_x_start - local_x_start);
  size = 2 * slice_in_nx;
  for (i=0; i < size; i++) {
    slice[ptr1 + i] = slice_in[i];
  }
}

void herm_destroy(herm_t * herm)
{
  fftw_free(herm->slice);
  fftw_free(herm->slice_out);
  fftw_free(herm->slice_in);
  free(herm->local_x_out_floor);
  free(herm->local_x_out_ceil);
  free(herm->local_x_in_floor);
  free(herm->local_x_in_ceil);
  free(herm->local_x_start);
  free(herm->local_nx);  
  free(herm);
}

dft_ptr_t max(dft_ptr_t i1, dft_ptr_t i2)
{
  return (i1 > i2) ? i1 : i2;
}

dft_ptr_t min(dft_ptr_t i1, dft_ptr_t i2)
{
  return (i1 < i2) ? i1 : i2;
}

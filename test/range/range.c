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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif /* HAVE_MPI */

#include "range.h"

typedef struct range_test_array_ range_test_array_t; 

struct range_test_array_
{
  double * data;
  range_t local_range;
};

int range_test_get_rank()
{
  int rank;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else /* DON'T HAVE_MPI */
  rank = 0;
#endif /* HAVE_MPI */
  return rank;
}

int range_test_get_nproc()
{
  int nproc;
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else /* DON'T HAVE_MPI */
  nproc = 1;
#endif /* HAVE_MPI */
  return nproc;
}

range_test_array_t range_test_setup(ptrdiff_t size)
{
  range_test_array_t a;
  int rank = range_test_get_rank();
  int nproc = range_test_get_nproc();
  ptrdiff_t block_size = size / nproc;
  a.local_range.start = rank * block_size;
  a.local_range.size = (rank == nproc - 1) ? size - a.local_range.start : block_size;
  a.data = (double *) malloc( sizeof(double) * a.local_range.size );
  for (ptrdiff_t i = 0; i < a.local_range.size; ++i)
  {
    a.data[i] = (double) (a.local_range.start + i);
  }
  return a;
}

int range_test_check(range_test_array_t a, range_t src_range, range_t dest_range, ptrdiff_t size)
{
  int rank = range_test_get_rank();
  int nproc = range_test_get_nproc();

  ptrdiff_t block_size = size / nproc;
  ptrdiff_t local_start = rank * block_size;
  for (ptrdiff_t i = 0; i < a.local_range.size; ++i)
  {
    double value = a.data[i];
    double expected_value;
    if (dest_range.start <= local_start + i && local_start + i < dest_range.start + dest_range.size)
    {
      expected_value = (double) (i + local_start - dest_range.start + src_range.start);
    }
    else
    {
      expected_value = (double) (local_start + i);
    }
    if (fabs(value - expected_value) > (value + 1.0) * 1e-6)
    {
      printf("error: rank=%d, global index=%ld, value=%f, expected=%f\n", rank, local_start + i, value, expected_value);
      return 1;
    }
  }
  return 0;
}

void range_test_destroy(range_test_array_t a)
{
  free(a.data);
}

int main(int argc, char **argv)
{
  if (argc != 2) {
    printf("Wrong number of arguments\n");
    return 1;
  }

  int size = atoi(argv[1]);

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif /* HAVE_MPI */

  range_t src_range;
  src_range.start = 0.6 * size;
  src_range.size = size - src_range.start - 1;
  range_t dest_range;
  dest_range.start = 1;
  dest_range.size = src_range.size;

  range_test_array_t test_array = range_test_setup(size);

  range_copy(test_array.data, src_range, test_array.data, dest_range, test_array.local_range);

  int err = range_test_check(test_array, src_range, dest_range, size);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif /* HAVE_MPI */

  return err;
}

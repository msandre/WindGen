#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif /* HAVE_MPI */

#include "range.h"

range_t range_create(ptrdiff_t start, ptrdiff_t size)
{
  range_t r;
  r.start = start;
  r.size = size;
  return r;
}

/**
 * range_intersect:
 * 
 * Compute the intersection of two ranges.
 */
range_t range_intersect(range_t r1, range_t r2)
{
  ptrdiff_t start1 = r1.start;
  ptrdiff_t start2 = r2.start;
  ptrdiff_t end1 = r1.start + r1.size;
  ptrdiff_t end2 = r2.start + r2.size;
  ptrdiff_t start = (start1 > start2) ? start1 : start2;
  ptrdiff_t end = (end1 < end2) ? end1 : end2;
  range_t intersection;
  intersection.start = start;
  intersection.size = (end > start) ? end - start : 0;
  return intersection;
}

range_t range_shift(range_t r, ptrdiff_t offset)
{
  range_t shifted;
  shifted.start = r.start + offset;
  shifted.size = r.size;
  return shifted;
}

/**
 * range_copy:
 * 
 * Copy range from source to destination.
 * 
 * src -- pointer to the local data field (a subset of the global data field)
 * 
 * src_range -- [src_range.start, src_range.start + src_range.size) is the
 *              sub-interval of the global data field to be copied from.
 * 
 * dest_range -- [dest_range.start, dest_range.start + dest_range.size) is the
 *               sub-interval of the global data field to be copied to.
 * 
 * local_range -- [local_range.start, local_range.start + local_range.size) is the
 *                sub-interval of the global data field pointed to by src.
 */
void range_copy(double * src, range_t src_range, double * dest, range_t dest_range, range_t local_range)
{
  assert(src_range.size == dest_range.size);

  int rank;
  int nproc;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else /* DON'T HAVE_MPI */
  rank = 0;
  nproc = 1;
#endif /* HAVE_MPI */

  ptrdiff_t * local_start = (ptrdiff_t *) malloc( sizeof(ptrdiff_t) * nproc );
  ptrdiff_t * local_size = (ptrdiff_t *) malloc( sizeof(ptrdiff_t) * nproc );

#ifdef HAVE_MPI
  if (sizeof(ptrdiff_t) == sizeof(int)) {
    MPI_Allgather(&local_range.start, 1, MPI_INT, local_start, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&local_range.size, 1, MPI_INT, local_size, 1, MPI_INT, MPI_COMM_WORLD);
  }
  else if (sizeof(ptrdiff_t) == sizeof(long)) {
    MPI_Allgather(&local_range.start, 1, MPI_LONG, local_start, 1, MPI_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&local_range.size, 1, MPI_LONG, local_size, 1, MPI_LONG, MPI_COMM_WORLD);
  }
  else if (sizeof(ptrdiff_t) == sizeof(long long)) {
    MPI_Allgather(&local_range.start, 1, MPI_LONG_LONG, local_start, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgather(&local_range.size, 1, MPI_LONG_LONG, local_size, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
  }
  else {
    abort();
  }
#else /* DON'T HAVE_MPI */
  local_start[0] = local_range.start;
  local_size[0] = local_range.size;
#endif /* HAVE_MPI */

  for (int send_proc = 0; send_proc < nproc; ++send_proc)
  {
    for (int recv_proc = 0; recv_proc < nproc; ++recv_proc)
    {
      range_t sender_range = range_create(local_start[send_proc], local_size[send_proc]);
      range_t sender_subrange = range_intersect(src_range, sender_range);
      range_t receiver_range = range_create(local_start[recv_proc], local_size[recv_proc]);
      range_t receiver_subrange = range_intersect(dest_range, receiver_range);
      range_t src_subrange = range_shift(
        range_intersect(
          range_shift(sender_subrange, -src_range.start),
          range_shift(receiver_subrange, -dest_range.start)
        ),
        src_range.start
      );
      if (src_subrange.size > 0)
      {
        range_t dest_subrange = range_shift(src_subrange, dest_range.start - src_range.start);
        ptrdiff_t src_offset = src_subrange.start-sender_range.start;
        ptrdiff_t dest_offset = dest_subrange.start-receiver_range.start;
        if (rank == send_proc && rank == recv_proc)
        {
          for (ptrdiff_t i = 0; i < src_subrange.size; ++i)
          {
            *(dest+dest_offset+i) = *(src+src_offset+i);
          }
        }
#ifdef HAVE_MPI
        else if (rank == send_proc && rank != recv_proc)
        {
          MPI_Send(src+src_offset, src_subrange.size, MPI_DOUBLE, recv_proc, send_proc, MPI_COMM_WORLD);
        }
        else if (rank != send_proc && rank == recv_proc)
        {
          MPI_Recv(dest+dest_offset, dest_subrange.size, MPI_DOUBLE, send_proc, send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
#endif /* HAVE_MPI */
      }
    }
  }

  free(local_start);
  free(local_size);
}
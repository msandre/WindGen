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

#include "rv.h"

#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif /* HAVE_MPI */

int main(int argc, char **argv)
{
  if (argc != 2) {
    printf("Wrong number of arguments\n");
    return 1;
  }

  int i, nl, n, rank, nproc;
  double * x, * y, * xg, * yg;
  FILE * fp;
  nl = atoi(argv[1]);

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else /* DON'T HAVE_MPI */  
  rank = 0;
  nproc = 1;
#endif /* HAVE_MPI */

  n = nproc * nl;

  x = (double *) malloc( sizeof(double) * nl );
  for (i=0; i < nl; i++) x[i] = 0.;
  y = (double *) malloc( sizeof(double) * nl );
  for (i=0; i < nl; i++) y[i] = 0.;

  rv_init(rank);

  rv_normal(x, nl);
  rv_normal(y, nl);
#ifdef HAVE_MPI
  if (rank == 0) {
    xg = (double *) malloc( sizeof(double) * n );
    for (i=0; i < n; i++) xg[i] = 0.;
    yg = (double *) malloc( sizeof(double) * n );
    for (i=0; i < n; i++) yg[i] = 0.;
  }
  else {
    xg = NULL;
    yg = NULL;
  }

  MPI_Gather(x, nl, MPI_DOUBLE, xg, nl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(y, nl, MPI_DOUBLE, yg, nl, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else /* DON'T HAVE_MPI */  
  xg = x; x = NULL;
  yg = y; y = NULL;
#endif /* HAVE_MPI */

  if (rank == 0) {
    fp = fopen("rv.res", "w");
    for (i=0; i < n; i++)
      fprintf(fp, "%f %f\n", xg[i], yg[i]);
    fclose(fp);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif /* HAVE_MPI */

  free(x);
  free(y);
  free(xg);
  free(yg);
}

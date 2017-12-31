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

#ifdef HAVE_MPI
#include <mpi.h>
#endif /* HAVE_MPI */
#include "wfs.h"

void init_line_data(wfs_t * sim, double ** purec, double ** pvrec, double ** pwrec)
{
  dft_ptr_t local_nx = sim->dft_x->local_nx;
  double * urec, * vrec, * wrec;
  dft_ptr_t i;

  urec = (double *) malloc( sizeof(double) * local_nx );
  vrec = (double *) malloc( sizeof(double) * local_nx );
  wrec = (double *) malloc( sizeof(double) * local_nx );
  
  for (i=0; i < local_nx; i++) {
    urec[i] = 0.;
    vrec[i] = 0.;
    wrec[i] = 0.;
  }

  (* purec) = urec;
  (* pvrec) = vrec;
  (* pwrec) = wrec;
}

void set_line_data(wfs_t * sim, dft_ptr_t j, dft_ptr_t k, double * urec, double * vrec, double * wrec)
{
  dft_ptr_t ny       = sim->box->ny;
  dft_ptr_t nz       = sim->box->nz;
  dft_ptr_t local_nx = sim->dft_x->local_nx;
  double *  fx       = sim->dft_x->field;
  double *  fy       = sim->dft_y->field;
  double *  fz       = sim->dft_z->field;
  dft_ptr_t i, iptr, jptr, ptr;
  
  jptr = j * (nz + 2);

  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz + 2);
    ptr  = iptr + jptr + k;

    urec[i] = fx[ptr];
    vrec[i] = fy[ptr];
    wrec[i] = fz[ptr];
  }
}

void finalize_line_data(wfs_t * sim, double * urec, double * vrec, double * wrec, int N)
{
  dft_ptr_t nx       = sim->box->nx;
  dft_ptr_t local_nx = sim->dft_x->local_nx;
  double    dx       = sim->box->lx / (double) (nx - 1);
  FILE * fp;
  char file[256];
  int p, rank, nproc;
  dft_ptr_t i, ptr, * local_nx_proc;
  double x, * gurec, * gvrec, *gwrec;

  local_nx_proc = NULL;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  local_nx_proc = (dft_ptr_t *) malloc( sizeof(dft_ptr_t) * nproc );
  if (sizeof(dft_ptr_t) == sizeof(int)) {
    MPI_Allgather(&local_nx, 1, MPI_INT, local_nx_proc, 1, MPI_INT, MPI_COMM_WORLD);
  }
  else if (sizeof(dft_ptr_t) == sizeof(long)) {
    MPI_Allgather(&local_nx, 1, MPI_LONG, local_nx_proc, 1, MPI_LONG, MPI_COMM_WORLD);
  }
  else if (sizeof(dft_ptr_t) == sizeof(long long)) {
    MPI_Allgather(&local_nx, 1, MPI_LONG_LONG, local_nx_proc, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
  }
  else {
    free(local_nx_proc);
  }

  if (rank == 0) {
    gurec = (double *) malloc( sizeof(double) * nx );
    gvrec = (double *) malloc( sizeof(double) * nx );
    gwrec = (double *) malloc( sizeof(double) * nx );

    for (i=0; i < local_nx; i++) {
      gurec[i] = urec[i];
      gvrec[i] = vrec[i];
      gwrec[i] = wrec[i];
    }
  }
  else {
    gurec = NULL;
    gvrec = NULL;
    gwrec = NULL;
  }

  ptr = local_nx_proc[0];
  for (p=1; p < nproc; p++) {
    if (rank == p) {
      MPI_Send(urec, local_nx, MPI_DOUBLE, 0, p, MPI_COMM_WORLD);
      MPI_Send(vrec, local_nx, MPI_DOUBLE, 0, p, MPI_COMM_WORLD);
      MPI_Send(wrec, local_nx, MPI_DOUBLE, 0, p, MPI_COMM_WORLD);
    }
    else if (rank == 0) {
      MPI_Recv(&gurec[ptr], local_nx_proc[p], MPI_DOUBLE, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&gvrec[ptr], local_nx_proc[p], MPI_DOUBLE, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&gwrec[ptr], local_nx_proc[p], MPI_DOUBLE, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    ptr = ptr + local_nx_proc[p];
  }
#else /* DON'T HAVE_MPI */
  rank = 0;
  nproc = 1;
  gurec = urec; urec = NULL;
  gvrec = vrec; vrec = NULL;
  gwrec = wrec; wrec = NULL;
#endif /* HAVE_MPI */

  sprintf(file, "simiu_%d.res", N);

  if (rank == 0) {
    x = 0.;
    fp = fopen(file, "w");
    for (i=0; i < nx; i++) {
      fprintf(fp, "%f %f %f %f\n", x, gurec[i], gvrec[i], gwrec[i]);
      x = x + dx;
    }
    fclose(fp);
  }

  free(local_nx_proc);
  free(urec);
  free(vrec);
  free(wrec);
  free(gurec);
  free(gvrec);
  free(gwrec);
}

int main(int argc, char *argv[])
{
  wfs_t * sim;
  double * urec, * vrec, * wrec;
  int i, N;
  dft_ptr_t ny, nz;
  char err[WFS_ERROR_SIZE];

  if (argc != 3)
    {
      printf("Wrong Arguments.\n");
      return 1;
    }
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  wfs_init_mpi();
#endif /* HAVE_MPI */

  sim = wfs_init(argv[1], err);
  if (sim == NULL) {
    printf("%s\n", err);
    return 1;
  }

  N = atoi(argv[2]) - 1;
  ny = sim->box->ny;
  nz = sim->box->nz;

  init_line_data(sim, &urec, &vrec, &wrec);
  wfs_generate_wind(sim);
  set_line_data(sim, ny/2, nz/2, urec, vrec, wrec);
  finalize_line_data(sim, urec, vrec, wrec, N);
  wfs_destroy(sim);
  
#ifdef HAVE_MPI
  wfs_finalize_mpi();
  MPI_Finalize();
#endif /* HAVE_MPI */

  return 0;
}

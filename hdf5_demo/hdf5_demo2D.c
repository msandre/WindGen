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
#include <stdio.h>
#endif /* HAVE_MPI */
#include "wfs.h"

int main(int argc, char *argv[])
{
  wfs_t * sim;
  char err[WFS_ERROR_SIZE];
  int rank = 0;

  // input file name expected as argument
  if (argc != 2)
    {
      printf("Wrong Arguments.\n");
      return 1;
    }

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  wfs_init_mpi();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif /* HAVE_MPI */

  // initialize wind field simulation
  sim = wfs_init2D(argv[1], err);
  if (sim == NULL) {
    printf("%s\n", err);
    return 1;
  }

  // generate wind
  wfs_generate_wind2D(sim);

  /* write data to HDF5 file
   *
   * This file can be viewed in ParaView if the built-in "reader" VisItPixieReader
   * is selected. This is useful for visualization of the velocity fluctuations.
   */
  hio_write2D("windgen2D.h5", sim);

  // destroy simulation
  wfs_destroy2D(sim);

  if(rank == 0) {
    printf("Finished!\n", rank);
  }

#ifdef HAVE_MPI
  wfs_finalize_mpi();
  MPI_Finalize();
#endif /* HAVE_MPI */

  return 0;
}

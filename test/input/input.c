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

#include "wfs.h"

int main(int argc, char *argv[])
{
  wfs_t * sim;
  char err[WFS_ERROR_SIZE];

  if (argc != 2)
    {
      printf("Wrong Arguments. Please give valid filename\n");
      return 1;
    }

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif /* HAVE_MPI */

  sim = wfs_init(argv[1], err);

  if (sim == NULL) {
    printf("%s\n", err);
    return 1;
  }

  printf("lx = %f\n", sim->box->lx);
  printf("ly = %f\n", sim->box->ly);
  printf("lz = %f\n", sim->box->lz);
  printf("nx = %ld\n", sim->box->nx);
  printf("ny = %ld\n", sim->box->ny);
  printf("nz = %ld\n", sim->box->nz);
  printf("umean = %f\n", sim->umean);
  printf("height = %f\n", sim->height);
  printf("roughness = %f\n", sim->roughness);
  printf("log roughness = %f\n", sim->log_roughness);
  printf("spectra = %s\n", sim->spectra);
  printf("conv = %d\n", sim->conv);
  printf("utau = %f\n", sim->utau);
  printf("gamma = %f\n", sim->gamma);

  wfs_destroy(sim);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif /* HAVE_MPI */

  return 0;
}

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

#include "conv.h"
#include <math.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  double a[3][3];
  double d[3];
  double v[3][3];

  if (argc != 10) {
    printf("Wrong number of arguments\n");
    return 1;
  }
  
  a[0][0] = atof(argv[1]);
  a[0][1] = atof(argv[2]);
  a[0][2] = atof(argv[3]);
  a[1][0] = atof(argv[4]);
  a[1][1] = atof(argv[5]);
  a[1][2] = atof(argv[6]);
  a[2][0] = atof(argv[7]);
  a[2][1] = atof(argv[8]);
  a[2][2] = atof(argv[9]);
  
  conv_eigen_value(a, d);
  conv_eigen_vector(a, v, d);

  printf("Wfs eigenvalues: (%17.5f, %17.5f, %17.5f)\n", d[0], d[1], d[2]);
  printf("Wfs eigenvector: (%17.5f, %17.5f, %17.5f)\n", v[0][0], v[1][0], v[2][0]);
  printf("Wfs eigenvector: (%17.5f, %17.5f, %17.5f)\n", v[0][1], v[1][1], v[2][1]);
  printf("Wfs eigenvector: (%17.5f, %17.5f, %17.5f)\n", v[0][2], v[1][2], v[2][2]);

  return 0;
}

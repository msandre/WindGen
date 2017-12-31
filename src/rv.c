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
#include <math.h>
#include <time.h>

static double z, maxran;
static double v[98];

/**
 * rv_init:
 * Initialize random number generator. Adapted from
 * W. H. Press et al., Numerical Recipes in C, Cambridge
 * University Press, 1988.
 */
void rv_init(int rank)
{
  unsigned int i, k;
  int j, s;

  s = (rank + 1) * time(NULL);

  i = 2;
  do {
    k = i;
    i = i << 1;
  } while(i);

  maxran = (double) k;
  srand(s);
  for (j=1; j < 98; j++) s = rand();
  for (j=1; j < 98; j++) v[j] = rand();
  z = rand();
}

/**
 * rv_uniform:
 * Generate a uniform random number in [0,1). Adapted from
 * W. H. Press et al., Numerical Recipes in C, Cambridge
 * University Press, 1988.
 */
double rv_uniform()
{
  int j;
  j = 1 + 97.0 * z / maxran;
  z = v[j];
  v[j] = (double) rand();
  return z / maxran;
}

/**
 * rv_normal:
 * Generate an array of iid Gaussian random numbers 
 * with zero mean and variance = 1/2.
 */
void rv_normal(double * a, int size)
{
  double sqrt2 = 1.414213562373095;
  double twopi = 6.283185307179586;
  int i, istart;
  double x, y, tmp;

  istart = size % 2;

  if (istart == 1) { // odd size
    x = rv_uniform();
    y = rv_uniform();
    a[0] = cos(twopi * x)  *  sqrt(-2. * log(1. - y)) / sqrt2;
  }

  for (i=istart; i < size; i=i+2) {
    x = rv_uniform();
    y = rv_uniform();
    tmp = sqrt(-2. * log(1. - y)) / sqrt2;
    a[i] = cos(twopi * x) * tmp;
    a[i+1] = sin(twopi * x) * tmp;
  }
}

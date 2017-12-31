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

/*
 * Author   : michael.andre@tum.de
 * Created  : 2013-07-18
 */

#ifndef __CONV_H__
#define __CONV_H__

// must be odd
#define CONV_SIZE 81

#include "wfs.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef double (* conv_func_ptr_t) (wfs_t *, double, double, double);

  void   conv_init            (wfs_t * sim);

  void   conv_execute         (wfs_t * sim, double kx, double ky, double kz, double * pfx, double * pfy, double * pfz);

  double conv_integrate       (wfs_t * sim, conv_func_ptr_t fptr, double kx, double ky, double kz);

  void   conv_decompose       (double a[][3]);
  
  void   conv_eigen_value     (double a[][3], double d[]);

  void   conv_eigen_vector    (double a[][3], double v[][3], double d[]);

  double phi11                (wfs_t * sim, double kx, double ky, double kz);
  double phi21                (wfs_t * sim, double kx, double ky, double kz);
  double phi31                (wfs_t * sim, double kx, double ky, double kz);
  double phi12                (wfs_t * sim, double kx, double ky, double kz);
  double phi22                (wfs_t * sim, double kx, double ky, double kz);
  double phi32                (wfs_t * sim, double kx, double ky, double kz);
  double phi13                (wfs_t * sim, double kx, double ky, double kz);
  double phi23                (wfs_t * sim, double kx, double ky, double kz);
  double phi33                (wfs_t * sim, double kx, double ky, double kz);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CONV_H__ */

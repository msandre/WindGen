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
 * Created  : 2013-06-03
 */

#ifndef __WFS_H__
#define __WFS_H__

#include <stdio.h>
#include <stdlib.h>
#include "def.h"
#include "dft.h"
#include "range.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef struct _wfs wfs_t;
  typedef struct _box box_t;
  
  struct _wfs
  {
    /*< public >*/
    box_t * box;
    double  umean, height, roughness;
    double  log_roughness;
    double  utau;
    char    spectra[10];
    int     conv;
    dft_t * dft_x;
    dft_t * dft_y;
    dft_t * dft_z;

    /*< private >*/
    double gamma;
    double LL;
    double iso_spec_coef;
  };

  struct _box
  {
    /*< public >*/
    double lx, ly, lz;
    dft_ptr_t nx, ny, nz;
  };

  void    wfs_init_mpi                 (void);
  void    wfs_finalize_mpi             (void);

  wfs_t * wfs_init                     (const char * file, char * err);

  int     is_2d                        (const wfs_t * sim);

  void    wfs_generate_wind            (wfs_t * sim);

  void    wfs_destroy                  (wfs_t * sim);

  void    wfs_apply_spectrum           (wfs_t * sim);

  void    wfs_apply_spectrum3D         (wfs_t * sim);

  void    wfs_apply_spectrum2D         (wfs_t * sim);

  void    wfs_sheared_spectrum         (wfs_t * sim, double kx, double ky, double kz, double * pkz0, double * pkk0, double * paxz, double * payz, double * pazz);

  void    wfs_apply_symmetry           (wfs_t * sim, dft_t * dft);

  double  hyp2f1_series                (double a, double b, double c, double z);
  
  double  hyp2f1                       (double z);

  double  ndim_eddy_time               (wfs_t * sim, double kk);

  void    hio_write                    (const char * file, wfs_t * sim);

  void    hio_fcopy                    (double * in, double * out, dft_ptr_t nx, dft_ptr_t ny, dft_ptr_t nz);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __WFS_H__ */

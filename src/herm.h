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
 * Created  : 2013-06-12
 */

#ifndef __HERM_H__
#define __HERM_H__

#include <stdlib.h>
#include "def.h"
#include "dft.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef struct _herm herm_t;

  struct _herm
  {
    /*< public >*/
    void * slice;
    dft_ptr_t nx, ny;
    int rank, nproc;
    dft_ptr_t * local_x_start;
    dft_ptr_t * local_nx;

    /*< private >*/
    void * slice_in, * slice_out;
    dft_ptr_t * local_x_out_floor;
    dft_ptr_t * local_x_out_ceil;
    dft_ptr_t * local_x_in_floor;
    dft_ptr_t * local_x_in_ceil;
  };

  herm_t *  herm_init           (dft_ptr_t nx, dft_ptr_t ny, dft_ptr_t local_x_start, dft_ptr_t local_nx);

  herm_t *  herm_init2D         (dft_ptr_t nx, dft_ptr_t local_x_start, dft_ptr_t local_nx);

  void      herm_execute        (herm_t * herm);

  void      herm_execute2D      (herm_t * herm);

  void      herm_destroy        (herm_t * herm);

  dft_ptr_t max                (dft_ptr_t i1, dft_ptr_t i2); 
  dft_ptr_t min                (dft_ptr_t i1, dft_ptr_t i2); 

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __HERM_H__ */

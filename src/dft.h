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
 * Created  : 2013-06-06
 */

#ifndef __DFT_H__
#define __DFT_H__

#include <stdlib.h>
#ifdef HAVE_MPI
#include <fftw3-mpi.h>
#else /* DON'T HAVE_MPI */
#include <fftw3.h>
#endif /* HAVE_MPI */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef ptrdiff_t dft_ptr_t;

  typedef struct _dft dft_t;
  typedef fftw_complex complex_t;

  struct _dft
  {
    /*< public >*/
    void * field;
    dft_ptr_t local_nx, local_x_start;

    /*< private >*/
    fftw_plan plan;
  };

  void    dft_init_mpi         (void);
  void    dft_finalize_mpi     (void);


  dft_t * dft_init             (dft_ptr_t nx, dft_ptr_t ny, dft_ptr_t nz);
  
  dft_t * dft_init2D           (dft_ptr_t nx, dft_ptr_t nz);

  void    dft_execute          (dft_t * dft);

  void    dft_destroy          (dft_t * dft);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __DFT_H__ */

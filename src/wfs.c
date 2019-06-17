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
#include "input.h"
#include "rv.h"
#include "herm.h"
#include "conv.h"
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifdef HAVE_HDF5
#include "hdf5.h"
#endif /* HAVE_HDF5 */

void wfs_init_mpi(void)
{
  dft_init_mpi();
}

void wfs_finalize_mpi(void)
{
  dft_finalize_mpi();
}

/**
 * wfs_init:
 * Initialize a wind field simulation from an input file.
 * 
 * 2D wind fields are generated for the input parameter ry=0. In that case ly
 * determines the interval for integrating out the y-component of the spectral
 * tensor.
 */
wfs_t * wfs_init(const char * file, char * err)
{
  int i, r, pos, rank;
  dft_ptr_t N;
  wfs_t * sim;

  char * fbuf = read_file(file, '#');
  if (fbuf == NULL) WFS_ERROR(err, "Input Error: Cannot open file");

  char * wfs_buf = reg_match(fbuf, "WindFieldSimulation\\s*\\{.*\\}");
  if (wfs_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid WindFieldSimulation block");

  char * box_buf = reg_match(wfs_buf, "(\\{|\\s)Box\\s*\\{.*\\}");
  if (box_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid Box block");

  char * lx_buf = reg_match(box_buf, "(\\{|\\s)lx\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (lx_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid valid lx");

  char * ly_buf = reg_match(box_buf, "(\\{|\\s)ly\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (ly_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid ly");

  char * rx_buf = reg_match(box_buf, "(\\{|\\s)rx\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)");
  if (rx_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid rx");

  char * ry_buf = reg_match(box_buf, "(\\{|\\s)ry\\s*=\\s*[0-9]+(\\}|\\s)");
  if (ry_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid ry");

  char * lz_buf = reg_match(box_buf, "(\\{|\\s)lz\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (lz_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid lz");

  char * rz_buf = reg_match(box_buf, "(\\{|\\s)rz\\s*=\\s*[1-9]+[0-9]*(\\}|\\s)");
  if (rz_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid rz");

  char * umean_buf = reg_match(wfs_buf, "(\\{|\\s)umean\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (umean_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid umean");

  char * height_buf = reg_match(wfs_buf, "(\\{|\\s)height\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (height_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid height");

  char * roughness_buf = reg_match(wfs_buf, "(\\{|\\s)roughness\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (roughness_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid roughness");

  char * log_roughness_buf = reg_match(wfs_buf, "(\\{|\\s)log\\s*roughness\\s*=\\s*\\+?[0-9]+(\\.[0-9]*)?(\\}|\\s)");
  if (log_roughness_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid log roughness");

  char * spectra_buf = reg_match(wfs_buf, "(\\{|\\s)spectra\\s*=\\s*[A-Za-z]+(\\}|\\s)");
  if (spectra_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid spectra");

  char * conv_buf = reg_match(wfs_buf, "(\\{|\\s)conv\\s*=\\s*(0|1)(\\}|\\s)");
  if (conv_buf == NULL) WFS_ERROR(err, "Input Error: Cannot find valid conv");

  sim = (wfs_t *) malloc( sizeof(wfs_t) );
  sim->box = (box_t *) malloc( sizeof(box_t) );
  sim->box->lx = atof(strchr(lx_buf,'=') + 1);
  sim->box->ly = atof(strchr(ly_buf,'=') + 1);
  sim->box->lz = atof(strchr(lz_buf,'=') + 1);
  r = atoi(strchr(rx_buf,'=') + 1);
  sim->box->nx = 1;
  for (i=0; i<r; i++) sim->box->nx = 2 * sim->box->nx;
  r = atoi(strchr(ry_buf,'=') + 1);
  sim->box->ny = 1;
  for (i=0; i<r; i++) sim->box->ny = 2 * sim->box->ny;
  r = atoi(strchr(rz_buf,'=') + 1);
  sim->box->nz = 1;
  for (i=0; i<r; i++) sim->box->nz = 2 * sim->box->nz;
  N = sim->box->nx * sim->box->ny * sim->box->nz;
  if (is_2d(sim))
  {
    sim->dft_x = dft_init2D(sim->box->nx, sim->box->nz);
    sim->dft_y = dft_init2D(sim->box->nx, sim->box->nz);
    sim->dft_z = dft_init2D(sim->box->nx, sim->box->nz);
  }
  else
  {
    sim->dft_x = dft_init3D(sim->box->nx, sim->box->ny, sim->box->nz);
    sim->dft_y = dft_init3D(sim->box->nx, sim->box->ny, sim->box->nz);
    sim->dft_z = dft_init3D(sim->box->nx, sim->box->ny, sim->box->nz);
  }
  
  sim->umean = atof(strchr(umean_buf,'=') + 1);
  sim->height = atof(strchr(height_buf,'=') + 1);
  sim->roughness = atof(strchr(roughness_buf,'=') + 1);
  sim->log_roughness = atof(strchr(log_roughness_buf,'=') + 1);

  if (sim->height <= 0.) WFS_ERROR(err, "Input Error: Invalid height value");
  if (sim->roughness <= 0.) WFS_ERROR(err, "Input Error: Invalid roughness value");
  if (sim->log_roughness <= 0.) WFS_ERROR(err, "Input Error: Invalid log roughness value");
  sim->utau = 0.41 * sim->umean / log(sim->height / sim->roughness);

  pos = 0;
  while(spectra_buf[pos]) {
    spectra_buf[pos] = tolower(spectra_buf[pos]);
    pos = pos + 1;
  }
  if (strstr(spectra_buf, "kaimal") != NULL) {
    strcpy(sim->spectra, "kaimal");
    sim->gamma = 3.9;
    sim->LL = pow(0.59 * sim->height, 2);
    sim->iso_spec_coef = 3.2 * sim->utau * sim->utau / pow(sim->height, 2./3.) * pow(sim->LL, 17./6.) / 12.566370614359172;
    sim->conv = atoi(strchr(conv_buf,'=') + 1);
  }
  else if (strstr(spectra_buf, "simiu") != NULL) {
    strcpy(sim->spectra, "simiu");
    sim->gamma = 3.8;
    sim->LL = pow(0.79 * sim->height, 2);
    sim->iso_spec_coef = 2.8 * sim->utau * sim->utau / pow(sim->height, 2./3.) * pow(sim->LL, 17./6.) / 12.566370614359172;
    sim->conv = atoi(strchr(conv_buf,'=') + 1);
  }
  else if (strstr(spectra_buf, "isotropic") != NULL) {
    strcpy(sim->spectra, "isotropic");
    sim->gamma = 0.;
    sim->LL = pow(sim->height, 2);
    sim->iso_spec_coef = sim->umean * sim->umean / pow(sim->height, 2./3.) * pow(sim->LL, 17./6.) / 12.566370614359172;
    sim->conv = 0;
  }
  else {
    wfs_destroy(sim);
    WFS_ERROR(err, "Input Error: Unknown spectra");
  }
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
#else /* DON'T HAVE_MPI */
  rank = 0;
#endif /* HAVE_MPI */

  rv_init(rank);
  if (is_2d(sim))
  {
    conv_init2D(sim);
  }
  else
  {
    conv_init3D(sim);
  }
  

  free(fbuf);
  free(wfs_buf);
  free(box_buf);
  free(lx_buf);
  free(ly_buf);
  free(rx_buf);
  free(ry_buf);
  free(lz_buf);
  free(rz_buf);
  free(umean_buf);
  free(height_buf);
  free(roughness_buf);
  free(log_roughness_buf);
  free(spectra_buf);
  free(conv_buf);

  return sim;
}

int is_2d(const wfs_t * sim)
{
  return sim->box->ny == 1;
}

/**
 * wfs_generate_wind:
 * Simulate wind field fluctuations with desired second-order properties.
 */
void wfs_generate_wind(wfs_t * sim)
{
  dft_ptr_t nx            = sim->box->nx;
  dft_ptr_t ny            = sim->box->ny;
  dft_ptr_t nz            = sim->box->nz;
  dft_ptr_t local_x_start = sim->dft_x->local_x_start;
  dft_ptr_t local_nx      = sim->dft_x->local_nx;
  double *  fx            = sim->dft_x->field;
  double *  fy            = sim->dft_y->field;
  double *  fz            = sim->dft_z->field;
  dft_ptr_t size;

  // fill arrays with random complex numbers such that E[ZZ*] = 1
  size = 2 * local_nx * ny * (nz/2 + 1);
  rv_normal(fx, size);
  if (!is_2d(sim))
  {
    rv_normal(fy, size);
  }
  rv_normal(fz, size);

  // apply spectral tensor matrices
  wfs_apply_spectrum(sim);

  // apply Hermitian symmetry
  wfs_apply_symmetry(sim, sim->dft_x);
  wfs_apply_symmetry(sim, sim->dft_z);

  // inverse dft
  dft_execute(sim->dft_x);
  dft_execute(sim->dft_z);

  if (!is_2d(sim))
  {
    wfs_apply_symmetry(sim, sim->dft_y);
    dft_execute(sim->dft_y);
  }
}

void wfs_destroy(wfs_t * sim)
{
  free(sim->box);
  dft_destroy(sim->dft_x);
  dft_destroy(sim->dft_y);
  dft_destroy(sim->dft_z);
  free(sim);
}

/**
 * wfs_apply_spectrum:
 * Apply the isotropic/anisotropic transformation matrix.
 */
void wfs_apply_spectrum(wfs_t * sim)
{
  if (is_2d(sim))
  {
    wfs_apply_spectrum2D(sim);
  }
  else
  {
    wfs_apply_spectrum3D(sim);
  }
  
}

void wfs_apply_spectrum3D(wfs_t * sim)
{
  double       twopi         = 6.283185307179586;
  dft_ptr_t    nx            = sim->box->nx;
  dft_ptr_t    ny            = sim->box->ny;
  dft_ptr_t    nz            = sim->box->nz;
  dft_ptr_t    local_x_start = sim->dft_x->local_x_start;
  dft_ptr_t    local_nx      = sim->dft_x->local_nx;
  double       lx            = sim->box->lx + sim->box->lx / (double) (nx - 1);
  double       ly            = sim->box->ly + sim->box->ly / (double) (ny - 1);
  double       lz            = sim->box->lz + sim->box->lz / (double) (nz - 1);
  double       dkx           = twopi / lx;
  double       dky           = twopi / ly;
  double       dkz           = twopi / lz;
  double       dkk           = dkx * dky * dkz;
  complex_t *  fcx           = sim->dft_x->field;
  complex_t *  fcy           = sim->dft_y->field;
  complex_t *  fcz           = sim->dft_z->field;
  int          conv          = sim->conv;
  double cx, cy, cz;
  double axz, ayz, azz;
  double kx, ky, kz, kz0, kkx, kky, kk, kk0, tmp;
  dft_ptr_t i, j, k, l, iptr, jptr, ptr;
  dft_ptr_t i_pos_start, i_neg_start, i_pos_end, i_neg_end;

  // kx = 0
  if (local_x_start == 0) {
    i_pos_start = 1;
    i = 0;
    kx = 0;

    ky = 0;
    iptr = i * ny * (nz/2 + 1);
    // positive y-frequencies
    for (j=0; j < ny/2 + 1; j++) {
      kz = 0;
      jptr = j * (nz/2 + 1);
      for (k=0; k < nz/2 + 1; k++) {
	ptr = iptr + jptr + k;

	for (l=0; l < 2; l++) {
	  // zero streamwise mean values
	  fcx[ptr][l] = 0.0;
	  fcy[ptr][l] = 0.0;
	  fcz[ptr][l] = 0.0;
	}
	kz = kz + dkz;
      }
      ky = ky + dky;
    }

    ky = -dky * (double) (ny/2 - 1);
    // negative y-frequencies
    for (j=ny/2 + 1; j < ny; j++) {
      kz = 0;
      jptr = j * (nz/2 + 1);
      for (k=0; k < nz/2 + 1; k++) {
	ptr = iptr + jptr + k;
	
	for (l=0; l < 2; l++) {
	  // zero streamwise mean values
	  fcx[ptr][l] = 0.0;
	  fcy[ptr][l] = 0.0;
	  fcz[ptr][l] = 0.0;
	}
	kz = kz + dkz;
      }
      ky = ky + dky;
    }
  }
  else {
    i_pos_start = 0;
  }

  i_pos_end   = min(local_nx, nx/2 + 1 - local_x_start);
  i_neg_start = max(nx/2 + 1 - local_x_start, 0);
  i_neg_end   = local_nx;

  // positive x-frequencies
  kx = dkx * (double) (local_x_start + i_pos_start);
  for (i=i_pos_start; i < i_pos_end; i++) {
    kkx = kx * kx;
    ky = 0;
    iptr = i * ny * (nz/2 + 1);
    // positive y-frequencies
    for (j=0; j < ny/2 + 1; j++) {
      kky = ky * ky;
      kz = 0;
      jptr = j * (nz/2 + 1);
      for (k=0; k < nz/2 + 1; k++) {
	ptr = iptr + jptr + k;
	kk = kkx + kky + kz * kz;
	if (conv == 1 && kk < 9.0 / sim->LL) {
	  conv_execute3D(sim, kx, ky, kz, &fcx[ptr][0], &fcy[ptr][0], &fcz[ptr][0]);
	}
	else {
	  wfs_sheared_spectrum(sim, kx, ky, kz, &kz0, &kk0, &axz, &ayz, &azz);
	  tmp = sqrt(sim->iso_spec_coef * dkk / pow(1. + sim->LL * kk0, 17./6.));

	  for (l=0; l < 2; l++) {
	    // isotropic tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];

	    fcx[ptr][l] = tmp * ( kz0*cy - ky*cz);
	    fcy[ptr][l] = tmp * (-kz0*cx + kx*cz);
	    fcz[ptr][l] = tmp * ( ky *cx - kx*cy);

	    // sheared tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	  
	    fcx[ptr][l] = cx + axz * cz;
	    fcy[ptr][l] = cy + ayz * cz;
	    fcz[ptr][l] = azz * cz;
	  }
	}
	kz = kz + dkz;
      }
      ky = ky + dky;
    }

    ky = -dky * (double) (ny/2 - 1);
    // negative y-frequencies
    for (j=ny/2 + 1; j < ny; j++) {
      kky = ky * ky;
      kz = 0;
      jptr = j * (nz/2 + 1);
      for (k=0; k < nz/2 + 1; k++) {
	ptr = iptr + jptr + k;
	kk = kkx + kky + kz * kz;
	if (conv == 1 && kk < 9.0 / sim->LL) {
	  conv_execute3D(sim, kx, ky, kz, &fcx[ptr][0], &fcy[ptr][0], &fcz[ptr][0]);
	}
	else {
	  wfs_sheared_spectrum(sim, kx, ky, kz, &kz0, &kk0, &axz, &ayz, &azz);
	  tmp = sqrt(sim->iso_spec_coef * dkk / pow(1. + sim->LL * kk0, 17./6.));
	
	  for (l=0; l < 2; l++) {
	    // isotropic tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	    
	    fcx[ptr][l] = tmp * ( kz0*cy - ky*cz);
	    fcy[ptr][l] = tmp * (-kz0*cx + kx*cz);
	    fcz[ptr][l] = tmp * ( ky *cx - kx*cy);
	    
	    // sheared tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	    
	    fcx[ptr][l] = cx + axz * cz;
	    fcy[ptr][l] = cy + ayz * cz;
	    fcz[ptr][l] = azz * cz;
	  }
	}
	kz = kz + dkz;
      }
      ky = ky + dky;
    }
    kx = kx + dkx;
  }

  // negative x-frequencies
  kx = -dkx * (double) ( nx - (local_x_start + i_neg_start) );
  for (i=i_neg_start; i < i_neg_end; i++) {
    kkx = kx * kx;
    ky = 0;
    iptr = i * ny * (nz/2 + 1);
    // positive y-frequencies
    for (j=0; j < ny/2 + 1; j++) {
      kky = ky * ky;
      kz = 0;
      jptr = j * (nz/2 + 1);
      for (k=0; k < nz/2 + 1; k++) {
	ptr = iptr + jptr + k;
	kk = kkx + kky + kz * kz;
	if (conv == 1 && kk < 9.0 / sim->LL) {
	  conv_execute3D(sim, kx, ky, kz, &fcx[ptr][0], &fcy[ptr][0], &fcz[ptr][0]);
	}
	else {
	  wfs_sheared_spectrum(sim, kx, ky, kz, &kz0, &kk0, &axz, &ayz, &azz);
	  tmp = sqrt(sim->iso_spec_coef * dkk / pow(1. + sim->LL * kk0, 17./6.));

	  for (l=0; l < 2; l++) {
	    // isotropic tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	  
	    fcx[ptr][l] = tmp * ( kz0*cy - ky*cz);
	    fcy[ptr][l] = tmp * (-kz0*cx + kx*cz);
	    fcz[ptr][l] = tmp * ( ky *cx - kx*cy);
	  
	    // sheared tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	  
	    fcx[ptr][l] = cx + axz * cz;
	    fcy[ptr][l] = cy + ayz * cz;
	    fcz[ptr][l] = azz * cz;
	  }
	}
	kz = kz + dkz;
      }
      ky = ky + dky;
    }

    ky = -dky * (double) (ny/2 - 1);
    // negative y-frequencies
    for (j=ny/2 + 1; j < ny; j++) {
      kky = ky * ky;
      kz = 0;
      jptr = j * (nz/2 + 1);
      for (k=0; k < nz/2 + 1; k++) {
	ptr = iptr + jptr + k;
	kk = kkx + kky + kz * kz;
	if (conv == 1 && kk < 9.0 / sim->LL) {
	  conv_execute3D(sim, kx, ky, kz, &fcx[ptr][0], &fcy[ptr][0], &fcz[ptr][0]);
	}
	else {
	  wfs_sheared_spectrum(sim, kx, ky, kz, &kz0, &kk0, &axz, &ayz, &azz);
	  tmp = sqrt(sim->iso_spec_coef * dkk / pow(1. + sim->LL * kk0, 17./6.));
	  
	  for (l=0; l < 2; l++) {
	    // isotropic tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	  
	    fcx[ptr][l] = tmp * ( kz0*cy - ky*cz);
	    fcy[ptr][l] = tmp * (-kz0*cx + kx*cz);
	    fcz[ptr][l] = tmp * ( ky *cx - kx*cy);
	  
	    // sheared tensor
	    cx = fcx[ptr][l];
	    cy = fcy[ptr][l];
	    cz = fcz[ptr][l];
	  
	    fcx[ptr][l] = cx + axz * cz;
	    fcy[ptr][l] = cy + ayz * cz;
	    fcz[ptr][l] = azz * cz;
	  }
	}
	kz = kz + dkz;
      }
      ky = ky + dky;
    }
    kx = kx + dkx;
  }  
}

void wfs_apply_spectrum2D(wfs_t * sim)
{
  double       twopi         = 6.283185307179586;
  dft_ptr_t    nx            = sim->box->nx;
  dft_ptr_t    nz            = sim->box->nz;
  dft_ptr_t    local_x_start = sim->dft_x->local_x_start;
  dft_ptr_t    local_nx      = sim->dft_x->local_nx;
  double       lx            = sim->box->lx + sim->box->lx / (double) (nx - 1);
  double       lz            = sim->box->lz + sim->box->lz / (double) (nz - 1);
  double       dkx           = twopi / lx;
  double       dkz           = twopi / lz;
  double       dkk           = dkx * dkz;
  complex_t *  fcx           = sim->dft_x->field;
  complex_t *  fcz           = sim->dft_z->field;
  int          conv          = sim->conv;
  double cx, cy, cz;
  double axz, ayz, azz;
  double kx, ky, kz, kz0, kkx, kky, kk, kk0, tmp;
  dft_ptr_t i, j, k, l, iptr, jptr, ptr;
  dft_ptr_t i_pos_start, i_neg_start, i_pos_end, i_neg_end;

  // kx = 0
  if (local_x_start == 0) {
    i_pos_start = 1;
    i = 0;
    kx = 0;
    kz = 0;
    iptr = i * (nz/2 + 1);
    for (k=0; k < nz/2 + 1; k++) {
      ptr = iptr + k;
      for (l=0; l < 2; l++) {
        // zero streamwise mean values
        fcx[ptr][l] = 0.0;
        fcz[ptr][l] = 0.0;
      }
      kz = kz + dkz;
    }
  }
  else {
    i_pos_start = 0;
  }

  i_pos_end   = min(local_nx, nx/2 + 1 - local_x_start);
  i_neg_start = max(nx/2 + 1 - local_x_start, 0);
  i_neg_end   = local_nx;

  // positive x-frequencies
  kx = dkx * (double) (local_x_start + i_pos_start);
  for (i=i_pos_start; i < i_pos_end; i++) {
    kkx = kx * kx;
    iptr = i * (nz/2 + 1);
    kz = 0;
    for (k=0; k < nz/2 + 1; k++) {
      ptr = iptr + k;
      kk = kkx + kz * kz;
      conv_execute2D(sim, kx, kz, &fcx[ptr][0], &fcz[ptr][0]);
      kz = kz + dkz;
    }
    kx = kx + dkx;
  }

  // negative x-frequencies
  kx = -dkx * (double) ( nx - (local_x_start + i_neg_start) );
  for (i=i_neg_start; i < i_neg_end; i++) {
    kkx = kx * kx;
    iptr = i * (nz/2 + 1);
    kz = 0;
    for (k=0; k < nz/2 + 1; k++) {
      ptr = iptr + k;
      kk = kkx + kz * kz;
      conv_execute2D(sim, kx, kz, &fcx[ptr][0], &fcz[ptr][0]);
      kz = kz + dkz;
    }
    kx = kx + dkx;
  }
}

/**
 * wfs_sheared_spectrum:
 * Calculate the sheared tensor coefficients.
 */
void wfs_sheared_spectrum(wfs_t * sim, double kx, double ky, double kz, double * pkz0, double * pkk0, double * paxz, double * payz, double * pazz)
{
  double kkx  = kx * kx;
  double kky  = ky * ky;
  double kk   = kkx + kky + kz * kz;
  double beta, kz0, kk0, c1, c2;
  double top, bot, theta;

  beta = ndim_eddy_time(sim, kk);
  kz0  = kz + beta * kx;
  kk0  = kkx + kky + kz0 * kz0;

  top = beta * kx * sqrt(kkx + kky);
  bot = (kk0 - kz0*kx*beta);

  theta = atan2(top, bot);

  c1 = beta * kkx * ( kk0 - 2.*kz0*kz0 + beta*kx*kz0 ) / ( kk * (kkx + kky) );
  c2 = ky * kk0 * theta / pow(kkx + kky, 1.5);

  (* pkz0) = kz0;
  (* pkk0) = kk0;
  (* paxz) = c1 - ky * c2 / kx;
  (* payz) = ky * c1 / kx + c2;
  (* pazz) = kk0 / kk;
}

/**
 * wfs_apply_symmetry:
 * Make the complex fields have Hermitian symmetry.
 */
void wfs_apply_symmetry(wfs_t * sim, dft_t * dft)
{
  dft_ptr_t    nx            = sim->box->nx;
  dft_ptr_t    ny            = sim->box->ny;
  dft_ptr_t    nz            = sim->box->nz;
  dft_ptr_t    local_x_start = sim->dft_x->local_x_start;
  dft_ptr_t    local_nx      = sim->dft_x->local_nx;  
  complex_t *  fc            = dft->field;
  dft_ptr_t i, j, k, l, iptr, jptr, ptr;
  complex_t * slice;
  herm_t * herm;

  herm = herm_init(nx, ny, local_x_start, local_nx);

  slice = herm->slice;

  // set Hermitian planes k = 0
  k = 0;
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        slice[i*ny + j][l] = fc[ptr][l];
      }
    }
  }

  // apply Hermitian symmetry k = 0
  herm_execute(herm);

  // copy Hermitian planes k = 0
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        fc[ptr][l] = slice[i*ny + j][l];
      }
    }
  }

  // set Hermitian planes k = nz/2
  k = nz/2;
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        slice[i*ny + j][l] = fc[ptr][l];
      }
    }
  }

  // apply Hermitian symmetry k = nz/2
  herm_execute(herm);

  // copy Hermitian planes k = nz/2
  for (i=0; i < local_nx; i++) {
    iptr = i * ny * (nz/2 + 1);
    for (j=0; j < ny; j++) {
      jptr = j * (nz/2 + 1);
      ptr = iptr + jptr + k;

      for (l=0; l < 2; l++) {
        fc[ptr][l] = slice[i*ny + j][l];
      }
    }
  }

  herm_destroy(herm);
}

/**
 * wfs_apply_symmetry2D:
 * Make the complex fields have Hermitian symmetry.
 */
void wfs_apply_symmetry2D(wfs_t * sim)
{
  dft_ptr_t    nx            = sim->box->nx;
  dft_ptr_t    ny            = sim->box->ny;
  dft_ptr_t    nz            = sim->box->nz;
  dft_ptr_t    local_x_start = sim->dft_x->local_x_start;
  dft_ptr_t    local_nx      = sim->dft_x->local_nx;  
  complex_t *  fcx           = sim->dft_x->field;
  complex_t *  fcz           = sim->dft_z->field;
  dft_ptr_t i, j, k, l, iptr, jptr, ptr;
  complex_t * slicex, * slicey, * slicez;
  herm_t * hermx, * hermy, * hermz;

  hermx = herm_init(nx, ny, local_x_start, local_nx);
  hermz = herm_init(nx, ny, local_x_start, local_nx);

  slicex = hermx->slice;
  slicez = hermz->slice;

  // set Hermitian planes k = 0
  k = 0;
  for (i=0; i < local_nx; i++) {
    ptr = i * (nz/2 + 1) + k;
    for (l=0; l < 2; l++) {
      slicex[i][l] = fcx[ptr][l];
      slicez[i][l] = fcz[ptr][l];
    }
  }

  // apply Hermitian symmetry k = 0
  herm_execute(hermx);
  herm_execute(hermz);

  // copy Hermitian planes k = 0
  for (i=0; i < local_nx; i++) {
    iptr = i * (nz/2 + 1) + k;
    for (l=0; l < 2; l++) {
      fcx[ptr][l] = slicex[i][l];
      fcz[ptr][l] = slicez[i][l];
    }
  }

  // set Hermitian planes k = nz/2
  k = nz/2;
  for (i=0; i < local_nx; i++) {
    iptr = i * (nz/2 + 1) + k;
    for (l=0; l < 2; l++) {
      slicex[i][l] = fcx[ptr][l];
      slicez[i][l] = fcz[ptr][l];
    }
  }

  // apply Hermitian symmetry k = nz/2
  herm_execute(hermx);
  herm_execute(hermz);

  // copy Hermitian planes k = nz/2
  for (i=0; i < local_nx; i++) {
    iptr = i * (nz/2 + 1) + k;
    for (l=0; l < 2; l++) {
      fcx[ptr][l] = slicex[i][l];
      fcz[ptr][l] = slicez[i][l];
    }
  }

  herm_destroy(hermx);
  herm_destroy(hermz);
}

/**
 * hyp2f1_series:
 * Approximate the Gauss hypergeometric series for |z| <= 0.5.
 */
double hyp2f1_series(double a, double b, double c, double z)
{
  int i;
  double val, coef, zz;

  zz = z;
  val = 1.;
  coef = a * b / c;
  for (i=1; i < 5; i++) {
    val = val + coef * zz;
    zz = zz * z;
    coef = coef * (a + (double) i) * (b + (double) i) / (c + (double) i) / (1. + (double) i);
  }

  val = val + coef * zz;

  return val;
}

/**
 * hyp2f1:
 * Compute the hypergeometric function for z < 0.
 */
double hyp2f1(double z)
{
  static double  a     = 0.3333333333333333;
  static double  b     = 2.8333333333333335;
  static double  c     = 1.3333333333333333;
  static double  gc1   = 0.68834394261431353;
  static double  gc2   =-0.13333333333333333;
  double w, a1, b1, c1, a2, b2, c2;

  if (z < -1.) {
    w = 1. / (1. - z);
    a1 = a;
    b1 = c - b;
    c1 = a - b + 1.;
    a2 = b;
    b2 = c - a;
    c2 = b - a + 1.;
    return pow(1. - z, -a) * gc1 * hyp2f1_series(a1, b1, c1, w)
         + pow(1. - z, -b) * gc2 * hyp2f1_series(a2, b2, c2, w);
  }
  else {
    w = z / (z - 1.);
    a1 = a;
    b1 = c - b;
    c1 = c;
    return pow(1. - z, -a) * hyp2f1_series(a1, b1, c1, w);
  }
}

/**
 * ndim_eddy_time:
 * Compute non-dimensional eddy lifetime.
 */
double ndim_eddy_time(wfs_t * sim, double kk)
{
  double kkLL = kk * sim->LL;
  return sim->gamma * pow(kkLL, -0.333333333) / sqrt( hyp2f1(-1. / kkLL) );
}

/**
 * hio_write:
 * Write simulation data in HDF5.
 */
void hio_write(const char * file, wfs_t * sim)
{
#ifdef HAVE_HDF5
  double *     fx              = sim->dft_x->field;
  double *     fy              = sim->dft_y->field;
  double *     fz              = sim->dft_z->field;
  dft_ptr_t    local_nx        = sim->dft_x->local_nx;
  dft_ptr_t    local_x_start   = sim->dft_x->local_x_start;
  dft_ptr_t    nx              = sim->box->nx;
  dft_ptr_t    ny              = sim->box->ny;
  dft_ptr_t    nz              = sim->box->nz;
  hid_t file_id, dset_id, dspace_id, dtype_id;
  hsize_t dim[3];
  int ndim;
  herr_t status;
  double * dset;

  ndim = (is_2d(sim)) ? 2 : 3;

#ifdef HAVE_MPI
  hsize_t count[3];
  hsize_t offset[3];
  hid_t fspace_id, plist_id;
  int rank, nproc;

  dset = (double *) malloc( sizeof(double) * local_nx * ny * nz );

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  status = H5Pclose(plist_id);

  // lx
  dim[0] = 1;
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  dspace_id = H5Screate_simple(1, dim, NULL);
  dset_id = H5Dcreate2(file_id, "/lx", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->box->lx;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);

  // ly
  dset_id = H5Dcreate2(file_id, "/ly", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->box->ly;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);

  // lz
  dset_id = H5Dcreate2(file_id, "/lz", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->box->lz;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);

  // z
  dset_id = H5Dcreate2(file_id, "/z", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->height;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);

  // umean
  dset_id = H5Dcreate2(file_id, "/umean", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->umean;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);
  
  // z0
  dset_id = H5Dcreate2(file_id, "/z0", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->roughness;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);

  // log_z0
  dset_id = H5Dcreate2(file_id, "/log_z0", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->log_roughness;
  if (rank == 0)
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, dset);
  status = H5Dclose(dset_id);

  // spectra
  dtype_id = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(dtype_id, strlen(sim->spectra) + 1);
  dset_id = H5Dcreate2(file_id, "/spectra", dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == 0)
    status = H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, plist_id, sim->spectra);
  status = H5Dclose(dset_id);
  status = H5Tclose(dtype_id);
  status = H5Sclose(dspace_id);
  status = H5Pclose(plist_id);

  dim[0] = nx;
  if (is_2d(sim))
  {
    dim[1] = nz;
  }
  else
  {
    dim[1] = ny;
    dim[2] = nz;
  }
  
  count[0] = local_nx;
  count[1] = dim[1];
  count[2] = dim[2];
  offset[0] = local_x_start;
  offset[1] = 0;
  offset[2] = 0;

  // u
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  fspace_id = H5Screate_simple(ndim, dim, NULL);
  dset_id = H5Dcreate2(file_id, "/u", H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sclose(fspace_id);
  dspace_id = H5Screate_simple(ndim, count, NULL);
  hio_fcopy(fx, dset, local_nx, ny, nz);
  fspace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, fspace_id, plist_id, dset);
  status = H5Sclose(fspace_id);
  status = H5Sclose(dspace_id);
  status = H5Dclose(dset_id);

  // v
  if (!is_2d(sim))
  {
    fspace_id = H5Screate_simple(ndim, dim, NULL);
    dset_id = H5Dcreate2(file_id, "/v", H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Sclose(fspace_id);
    dspace_id = H5Screate_simple(ndim, count, NULL);
    hio_fcopy(fy, dset, local_nx, ny, nz);
    fspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, fspace_id, plist_id, dset);
    status = H5Sclose(fspace_id);
    status = H5Sclose(dspace_id);
    status = H5Dclose(dset_id);
  }

  // w
  fspace_id = H5Screate_simple(ndim, dim, NULL);
  dset_id = H5Dcreate2(file_id, "/w", H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sclose(fspace_id);
  dspace_id = H5Screate_simple(ndim, count, NULL);
  hio_fcopy(fz, dset, local_nx, ny, nz);
  fspace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, fspace_id, plist_id, dset);
  status = H5Sclose(fspace_id);
  status = H5Sclose(dspace_id);
  status = H5Dclose(dset_id);
  status = H5Pclose(plist_id);
  status = H5Fclose(file_id);
  free(dset);

#else /* DON'T HAVE_MPI */

  dset = (double *) malloc( sizeof(double) * local_nx * ny * nz );

  file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // lx
  dim[0] = 1;
  dspace_id = H5Screate_simple(1, dim, NULL);
  dset_id = H5Dcreate2(file_id, "/lx", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->box->lx;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // ly
  dset_id = H5Dcreate2(file_id, "/ly", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->box->ly;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // lz
  dset_id = H5Dcreate2(file_id, "/lz", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->box->lz;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // z
  dset_id = H5Dcreate2(file_id, "/z", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->height;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // umean
  dset_id = H5Dcreate2(file_id, "/umean", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->umean;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // z0
  dset_id = H5Dcreate2(file_id, "/z0", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->roughness;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // log_z0
  dset_id = H5Dcreate2(file_id, "/log_z0", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dset[0] = sim->log_roughness;
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // spectra
  dtype_id = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(dtype_id, strlen(sim->spectra) + 1);
  dset_id = H5Dcreate2(file_id, "/spectra", dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, sim->spectra);
  status = H5Dclose(dset_id);
  status = H5Tclose(dtype_id);
  status = H5Sclose(dspace_id);

  dim[0] = nx;
  if (is_2d(sim))
  {
    dim[1] = nz;
  }
  else
  {
    dim[1] = ny;
    dim[2] = nz;
  }
  dspace_id = H5Screate_simple(ndim, dim, NULL);

  // u
  dset_id = H5Dcreate2(file_id, "/u", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hio_fcopy(fx, dset, nx, ny, nz);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);

  // v
  if (!is_2d(sim))
  {
    dset_id = H5Dcreate2(file_id, "/v", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hio_fcopy(fy, dset, nx, ny, nz);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
    status = H5Dclose(dset_id);
  }

  // w
  dset_id = H5Dcreate2(file_id, "/w", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hio_fcopy(fz, dset, nx, ny, nz);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset);
  status = H5Dclose(dset_id);
  status = H5Sclose(dspace_id);
  status = H5Fclose(file_id);
  free(dset);
#endif /* HAVE_MPI */
#endif /* HAVE_HDF5 */
}

/**
 * hio_fcopy:
 * Copy padded DFT field to non-padded field.
 */
void hio_fcopy(double * in, double * out, dft_ptr_t local_nx, dft_ptr_t ny, dft_ptr_t nz)
{
#ifdef HAVE_HDF5
  dft_ptr_t i, j, k, iptr1, iptr2, jptr1, jptr2;

  for (i=0; i < local_nx; i++) {
    iptr1 = i * ny * nz;
    iptr2 = i * ny * (nz + 2);
    for (j=0; j < ny; j++) {
      jptr1 = j * nz;
      jptr2 = j * (nz + 2);
      for (k=0; k < nz; k++)
	out[iptr1 + jptr1 + k] = in[iptr2 + jptr2 + k];
    }
  }
#endif /* HAVE_HDF5 */
}

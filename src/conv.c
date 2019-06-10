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

static double fourpi = 12.566370614359172;
static double sqrt3  = 1.7320508075688772;
static double dsy, dsz, ds;
static double sinc2y[CONV_SIZE];
static double sinc2z[CONV_SIZE];

/**
 * conv_init:
 * Initialize the convolution arrays.
 */
void conv_init(wfs_t * sim)
{
  double   lx    = sim->box->lx;
  double   ly    = sim->box->ly;
  double   lz    = sim->box->lz;
  int      n     = CONV_SIZE;
  int i;
  double x;

  dsy = fourpi / ((double) n * ly);
  dsz = fourpi / ((double) n * lz);
  ds = fourpi / (2 * lx) * 1.1076 * 1.1076 * dsy * dsz;
  
  // y-component
  for (i = -(n-1)/2; i <= (n-1)/2; i++) {
    if (i == 0) {
      sinc2y[i + (n-1)/2] = 1.;
    }
    else {
      x = (double) -i * dsy * ly / 2.;
      sinc2y[i + (n-1)/2] = sin(x) * sin(x) / (x * x);
    }
  }

  // z-component
  for (i = -(n-1)/2; i <= (n-1)/2; i++) {
    if (i == 0) {
      sinc2z[i + (n-1)/2] = 1.;
    }
    else {
      x = (double) -i * dsz * lz / 2.;
      sinc2z[i + (n-1)/2] = sin(x) * sin(x) / (x * x);
    }
  }
}

/**
 * conv_init2D:
 * Initialize the convolution arrays.
 */
void conv_init2D(wfs_t * sim)
{
  double   lx    = sim->box->lx;
  double   lz    = sim->box->lz;
  int      n     = CONV_SIZE;
  int i;
  double x, sinx;

  dsz = fourpi / ((double) n * lz);
  ds = fourpi / (2 * lx) * 1.1076 * dsz;
  
  // z-component
  for (i = -(n-1)/2; i <= (n-1)/2; i++) {
    if (i == 0) {
      sinc2z[i + (n-1)/2] = 1.;
    }
    else {
      x = (double) -i * dsz * lz / 2.;
      sinc2z[i + (n-1)/2] = sin(x) * sin(x) / (x * x);
    }
  }
}

/**
 * conv_execute:
 * Substitute the dirac approximation with the more accurate convolution decomposition.
 */
void conv_execute(wfs_t * sim, double kx, double ky, double kz, double * pfx, double * pfy, double * pfz)
{
  double a[3][3];
  double cx, cy, cz;

  a[0][0] = conv_integrate(sim, &phi11, kx, ky, kz);
  a[1][0] = conv_integrate(sim, &phi21, kx, ky, kz);
  a[2][0] = conv_integrate(sim, &phi31, kx, ky, kz);
  a[0][1] = a[1][0];
  a[1][1] = conv_integrate(sim, &phi22, kx, ky, kz);
  a[2][1] = conv_integrate(sim, &phi32, kx, ky, kz);
  a[0][2] = a[2][0];
  a[1][2] = a[2][1];
  a[2][2] = conv_integrate(sim, &phi33, kx, ky, kz);
    
  conv_decompose(a);
  
  // real
  cx = pfx[0];
  cy = pfy[0];
  cz = pfz[0];
  pfx[0] = a[0][0] * cx + a[0][1] * cy + a[0][2] * cz;
  pfy[0] = a[1][0] * cx + a[1][1] * cy + a[1][2] * cz;
  pfz[0] = a[2][0] * cx + a[2][1] * cy + a[2][2] * cz;
    
  // imaginary
  cx = pfx[1];
  cy = pfy[1];
  cz = pfz[1];
  pfx[1] = a[0][0] * cx + a[0][1] * cy + a[0][2] * cz;
  pfy[1] = a[1][0] * cx + a[1][1] * cy + a[1][2] * cz;
  pfz[1] = a[2][0] * cx + a[2][1] * cy + a[2][2] * cz;
}

/**
 * conv_execute2D:
 * Substitute the dirac approximation with the more accurate convolution decomposition.
 */
void conv_execute2D(wfs_t * sim, double kx, double kz, double * pfx, double * pfz)
{
  double a[2][2];
  double cx, cz;

  a[0][0] = conv_integrate2D(sim, &phi11, kx, kz);
  a[1][0] = conv_integrate2D(sim, &phi31, kx, kz);
  a[0][1] = a[1][0];
  a[1][1] = conv_integrate2D(sim, &phi33, kx, kz);
    
  conv_decompose2D(a);
  
  // real
  cx = pfx[0];
  cz = pfz[0];
  pfx[0] = a[0][0] * cx + a[0][1] * cz;
  pfz[0] = a[1][0] * cx + a[1][1] * cz;
    
  // imaginary
  cx = pfx[1];
  cz = pfz[1];
  pfx[1] = a[0][0] * cx + a[0][1] * cz;
  pfz[1] = a[1][0] * cx + a[1][1] * cz;
}

/**
 * conv_integrate:
 * Compute approximate convolution of spectral tensor component.
 */
double conv_integrate(wfs_t * sim, conv_func_ptr_t fptr, double kx, double ky, double kz)
{
  int n = CONV_SIZE;
  int iy, iz;
  double sy, sz, val;

  val = 0.;

  for (iy = -(n-1)/2; iy <= (n-1)/2; iy++) {
    sy = ky + (double) iy * dsy;
    for (iz = -(n-1)/2; iz <= (n-1)/2; iz++) {
      sz = kz + (double) iz * dsz;
      
      val += (* fptr) (sim, kx, sy, sz) * sinc2y[iy + (n-1)/2] * sinc2z[iz + (n-1)/2];
    }
  }
  
  return val * ds;
}

/**
 * conv_integrate2D:
 * Compute approximate convolution of spectral tensor component.
 * ly is used as the length over which to perform to the y-component integration.
 * CONV_SIZE is the number of integration points.
 */
double conv_integrate2D(wfs_t * sim, conv_func_ptr_t fptr, double kx, double kz)
{
  int n = CONV_SIZE;
  int iy, iz;
  double sy, sz, val;
  double ly = sim->box->ly;
  val = 0.;

  for (iy = -(n-1)/2; iy <= (n-1)/2; iy++) {
    sy = ((double) iy * 2 * ly ) / ((double) n - 1);
    for (iz = -(n-1)/2; iz <= (n-1)/2; iz++) {
      sz = kz + (double) iz * dsz;
      
      val += (* fptr) (sim, kx, sy, sz) * sinc2z[iz + (n-1)/2];
    }
  }
  
  return val * ds;
}

/**
 * conv_decompose:
 * Decompose a symmetric 3x3 matrix neglecting negative eigenvalues.
 */
void conv_decompose(double a[][3])
{
  double d[3];
  double v[3][3];

  conv_eigen_value(a, d);
  conv_eigen_vector(a, v, d);

  d[0] = (d[0] < 0.0) ? 0.0 : sqrt(d[0]);
  d[1] = (d[1] < 0.0) ? 0.0 : sqrt(d[1]);
  d[2] = (d[2] < 0.0) ? 0.0 : sqrt(d[2]);
  
  a[0][0] = d[0] * v[0][0];
  a[1][0] = d[0] * v[1][0];
  a[2][0] = d[0] * v[2][0];
  a[0][1] = d[1] * v[0][1];
  a[1][1] = d[1] * v[1][1];
  a[2][1] = d[1] * v[2][1];
  a[0][2] = d[2] * v[0][2];
  a[1][2] = d[2] * v[1][2];
  a[2][2] = d[2] * v[2][2];
}

/**
 * conv_decompose2D:
 * Decompose a symmetric 2x2 matrix neglecting negative eigenvalues.
 */
void conv_decompose2D(double a[][2])
{
  double tr, det, sdet, denom;
  double mp[2][2];

  tr    = a[0][0] + a[1][1];
  det   = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  sdet  = (det < 0.0) ? 0.0 : sqrt(det);
  denom = sqrt(tr + 2. * sdet);

  a[0][0] = (a[0][0] + sdet) * (1. / denom);
  a[0][1] = a[0][1] * (1. / denom);
  a[1][0] = a[0][1];
  a[1][1] = (a[1][1] + sdet) * (1. / denom);
}

/**
 * conv_eigen_value:
 * Compute the eigenvalues of a symmetric 3x3 matrix.
 */
void conv_eigen_value(double a[][3], double d[])
{
  double m, q, p, t, b11, b22, b33, sqrtp, cost, sint;
  
  // m = trace(a) / 3
  m = (a[0][0] + a[1][1] + a[2][2]) / 3.0;
  
  // q = det(a - m*I) / 2
  b11 = a[0][0] - m;
  b22 = a[1][1] - m;
  b33 = a[2][2] - m;
  q = 0.5 * (b11*b22*b33 - b11*a[1][2]*a[1][2] - b22*a[0][2]*a[0][2] - b33*a[0][1]*a[0][1] + 2.0*a[0][1]*a[0][2]*a[1][2]);
  
  // p = (a - m*I)_{i,j}*(a - m*I)_{i,j} / 6
  p = (b11*b11 + b22*b22 + b33*b33 + 2.0*(a[0][1]*a[0][1] + a[0][2]*a[0][2] + a[1][2]*a[1][2])) / 6.0;

  t = atan2(sqrt(p*p*p - q*q), q) / 3.0;
  sqrtp = sqrt(p);
  cost  = cos(t);
  sint  = sin(t);

  // eigenvalues
  d[0] = m + 2.0 * sqrtp * cost;
  d[1] = m - sqrtp * (cost + sqrt3 * sint);
  d[2] = m - sqrtp * (cost - sqrt3 * sint);
}

/**
 * conv_eigen_vector:
 * Compute the eigenvectors of a symmetric 3x3 matrix and its eigenvalues.
 */
void conv_eigen_vector(double a[][3], double v[][3], double d[])
{
  double c[3][3];
  int idx[3];
  int jdx[3];
  int i, j, k, l, imx, jmx, tmp;
  double vmax, fac, rdet;

  for (l=0; l < 3; l++) {
    // initialize
    for (i=0; i < 3; i++) {
      for (j=0; j < 3; j++) {
	c[i][j] = a[i][j];
      }
      c[i][i] -= d[l];
      idx[i] = i;
      jdx[i] = i;
    }

    for (k=0; k < 2; k++) {
      // move max |c_{i,j}| to idx[k], jdx[k]
      vmax = 0.0;
      for (i=k; i < 3; i++) {
	for (j=k; j < 3; j++) {
	  if (vmax < fabs(c[idx[i]][jdx[j]])) {
	    imx = i;
	    jmx = j;
	    vmax = fabs(c[idx[i]][jdx[j]]);
	  }
	}
      }
      tmp      = idx[k];
      idx[k]   = idx[imx];
      idx[imx] = tmp;
      tmp      = jdx[k];
      jdx[k]   = jdx[jmx];
      jdx[jmx] = tmp;

      // normalize row idx[k]
      fac = 1.0 / c[idx[k]][jdx[k]];
      c[idx[k]][jdx[k]] = 1.0;
      for (j=k+1; j < 3; j++)
	c[idx[k]][jdx[j]] *= fac;

      // eliminate non-zero entries in lower diagonal of column jdx[k]
      for (i=k+1; i < 3; i++) {
	fac = c[idx[i]][jdx[k]];
	c[idx[i]][jdx[k]] = 0.0;
	for (j=k+1; j < 3; j++)
	  c[idx[i]][jdx[j]] -= fac * c[idx[k]][jdx[j]];
      }
    }

    // solve for the lth eigenvector
    rdet = 1.0 / (c[idx[0]][jdx[1]]*c[idx[1]][jdx[2]] - c[idx[1]][jdx[1]]*c[idx[0]][jdx[2]]);
    v[jdx[0]][l] =-1.0;
    v[jdx[1]][l] = rdet * c[idx[1]][jdx[2]];
    v[jdx[2]][l] =-rdet * c[idx[1]][jdx[1]];
    
    // normalize
    fac = 1.0 / sqrt(v[0][l]*v[0][l] + v[1][l]*v[1][l] + v[2][l]*v[2][l]);
    for (j=0; j < 3; j++)
      v[j][l] *= fac;
  }
}

/**
 * phi11:
 * Component of the spectral tensor.
 */
double phi11(wfs_t * sim, double kx, double ky, double kz)
{
  double kkx    = kx * kx;
  double kky    = ky * ky;
  double kkz    = kz * kz;
  double kk     = kkx + kky + kkz;
  double beta   = sim->gamma * pow(kk * sim->LL, -0.333333) / sqrt(hyp2f1(-1. / (kk * sim->LL)));
  double kz0    = kz + beta * kx;
  double kk0    = kkx + kky + kz0 * kz0;
  double top    = beta * kx * sqrt(kkx + kky);
  double bot    = kk0 - kz0 * kx * beta;
  double c1     = beta * kkx * ( kk0 - 2. * kz0 * kz0 + beta * kx * kz0 ) / ( kk * (kkx + kky) );
  double c2     = ky * kk0 * atan2(top,bot) / pow(kkx + kky, 1.5);
  double zeta1  = c1 - ky * c2 / kx;
  double zeta2  = ky * c1 / kx + c2;
  double coef   = sim->iso_spec_coef / pow(1. + sim->LL * kk0, 17./6.);
  
  return coef * (kky * (zeta1 * zeta1 + 1.) + pow(kz0 - kx * zeta1, 2.));
}

/**
 * phi21:
 * Component of the spectral tensor.
 */
double phi21(wfs_t * sim, double kx, double ky, double kz)
{
  double kkx    = kx * kx;
  double kky    = ky * ky;
  double kkz    = kz * kz;
  double kk     = kkx + kky + kkz;
  double beta   = sim->gamma * pow(kk * sim->LL, -0.333333) / sqrt(hyp2f1(-1. / (kk * sim->LL)));
  double kz0    = kz + beta * kx;
  double kk0    = kkx + kky + kz0 * kz0;
  double top    = beta * kx * sqrt(kkx + kky);
  double bot    = kk0 - kz0 * kx * beta;
  double c1     = beta * kkx * ( kk0 - 2. * kz0 * kz0 + beta * kx * kz0 ) / ( kk * (kkx + kky) );
  double c2     = ky * kk0 * atan2(top,bot) / pow(kkx + kky, 1.5);
  double zeta1  = c1 - ky * c2 / kx;
  double zeta2  = ky * c1 / kx + c2;
  double coef   = sim->iso_spec_coef / pow(1. + sim->LL * kk0, 17./6.);
  
  return coef * (ky * zeta1 * (-kz0 + ky * zeta2) - kx * zeta2 * (kz0 - kx * zeta1) - kx * ky);
}

/**
 * phi31:
 * Component of the spectral tensor.
 */
double phi31(wfs_t * sim, double kx, double ky, double kz)
{
  double kkx    = kx * kx;
  double kky    = ky * ky;
  double kkz    = kz * kz;
  double kk     = kkx + kky + kkz;
  double beta   = sim->gamma * pow(kk * sim->LL, -0.333333) / sqrt(hyp2f1(-1. / (kk * sim->LL)));
  double kz0    = kz + beta * kx;
  double kk0    = kkx + kky + kz0 * kz0;
  double top    = beta * kx * sqrt(kkx + kky);
  double bot    = kk0 - kz0 * kx * beta;
  double c1     = beta * kkx * ( kk0 - 2. * kz0 * kz0 + beta * kx * kz0 ) / ( kk * (kkx + kky) );
  double c2     = ky * kk0 * atan2(top,bot) / pow(kkx + kky, 1.5);
  double zeta1  = c1 - ky * c2 / kx;
  double zeta2  = ky * c1 / kx + c2;
  double coef   = sim->iso_spec_coef / pow(1. + sim->LL * kk0, 17./6.);
  
  return coef * (kky * zeta1 - kx * (kz0 - kx * zeta1)) * kk0 / kk;
}

/**
 * phi12:
 * Component of the spectral tensor.
 */
double phi12(wfs_t * sim, double kx, double ky, double kz)
{
  return phi21(sim, kx, ky, kz);
}

/**
 * phi22:
 * Component of the spectral tensor.
 */
double phi22(wfs_t * sim, double kx, double ky, double kz)
{
  double kkx    = kx * kx;
  double kky    = ky * ky;
  double kkz    = kz * kz;
  double kk     = kkx + kky + kkz;
  double beta   = sim->gamma * pow(kk * sim->LL, -0.333333) / sqrt(hyp2f1(-1. / (kk * sim->LL)));
  double kz0    = kz + beta * kx;
  double kk0    = kkx + kky + kz0 * kz0;
  double top    = beta * kx * sqrt(kkx + kky);
  double bot    = kk0 - kz0 * kx * beta;
  double c1     = beta * kkx * ( kk0 - 2. * kz0 * kz0 + beta * kx * kz0 ) / ( kk * (kkx + kky) );
  double c2     = ky * kk0 * atan2(top,bot) / pow(kkx + kky, 1.5);
  double zeta1  = c1 - ky * c2 / kx;
  double zeta2  = ky * c1 / kx + c2;
  double coef   = sim->iso_spec_coef / pow(1. + sim->LL * kk0, 17./6.);
  
  return coef * (pow(-kz0 + ky * zeta2, 2.) + kkx * (zeta2 * zeta2 + 1.));
}

/**
 * phi32:
 * Component of the spectral tensor.
 */
double phi32(wfs_t * sim, double kx, double ky, double kz)
{
  double kkx    = kx * kx;
  double kky    = ky * ky;
  double kkz    = kz * kz;
  double kk     = kkx + kky + kkz;
  double beta   = sim->gamma * pow(kk * sim->LL, -0.333333) / sqrt(hyp2f1(-1. / (kk * sim->LL)));
  double kz0    = kz + beta * kx;
  double kk0    = kkx + kky + kz0 * kz0;
  double top    = beta * kx * sqrt(kkx + kky);
  double bot    = kk0 - kz0 * kx * beta;
  double c1     = beta * kkx * ( kk0 - 2. * kz0 * kz0 + beta * kx * kz0 ) / ( kk * (kkx + kky) );
  double c2     = ky * kk0 * atan2(top,bot) / pow(kkx + kky, 1.5);
  double zeta1  = c1 - ky * c2 / kx;
  double zeta2  = ky * c1 / kx + c2;
  double coef   = sim->iso_spec_coef / pow(1. + sim->LL * kk0, 17./6.);
  
  return coef * (-kz0 * ky + zeta2 * (kkx + kky)) * kk0 / kk;
}

/**
 * phi13:
 * Component of the spectral tensor.
 */
double phi13(wfs_t * sim, double kx, double ky, double kz)
{
  return phi31(sim, kx, ky, kz);
}

/**
 * phi23:
 * Component of the spectral tensor.
 */
double phi23(wfs_t * sim, double kx, double ky, double kz)
{
  return phi32(sim, kx, ky, kz);
}

/**
 * phi33:
 * Component of the spectral tensor.
 */
double phi33(wfs_t * sim, double kx, double ky, double kz)
{
  double kkx    = kx * kx;
  double kky    = ky * ky;
  double kkz    = kz * kz;
  double kk     = kkx + kky + kkz;
  double beta   = sim->gamma * pow(kk * sim->LL, -0.333333) / sqrt(hyp2f1(-1. / (kk * sim->LL)));
  double kz0    = kz + beta * kx;
  double kk0    = kkx + kky + kz0 * kz0;
  double top    = beta * kx * sqrt(kkx + kky);
  double bot    = kk0 - kz0 * kx * beta;
  double c1     = beta * kkx * ( kk0 - 2. * kz0 * kz0 + beta * kx * kz0 ) / ( kk * (kkx + kky) );
  double c2     = ky * kk0 * atan2(top,bot) / pow(kkx + kky, 1.5);
  double zeta1  = c1 - ky * c2 / kx;
  double zeta2  = ky * c1 / kx + c2;
  double coef   = sim->iso_spec_coef / pow(1. + sim->LL * kk0, 17./6.);
  
  return coef * (kkx + kky) * kk0 * kk0 / kk / kk;
}

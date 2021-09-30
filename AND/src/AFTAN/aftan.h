/*
This file is part of ANCC.

AFTAN is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AFTAN is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef AFTAN_H
#define AFTAN_H 1

/* Finction prorotypes */


void aftanipg_(double *PIover4, int *n, float *seis, float *t0, float *dt,
    float *delta, float *vmin, float *vmax, float *tmin, float *tmax,
    float *tresh, float *ffact, float *perc, int *npoints,
    float *taperl, int *nfin, float *snr, float *fmatch,
    int *npred, float **pred,
    int *nprpv, float *prpvper, float *prpvvel,
    int *nfout1, float **arr1, int *nfout2, float *arr2,
    float *tamp, int *nrow, int *ncol, float **ampo, int *ierr);

void aftanpg_(double *PIover4, int *n, float *seis, float *t0, float *dt,
    float *delta, float *vmin, float *vmax, float *tmin, float *tmax,
    float *tresh, float *ffact, float *perc, int *npoints,
    float *taperl, int *nfin, float *snr,
    int *nprpv, float *prpvper, float *prpvvel,
    int *nfout1, float **arr1, int *nfout2, float *arr2,
    float *tamp, int *nrow, int *ncol, float **ampo, int *ierr);

void printres_(float *dt, float *delta, int *nfout1, float **arr1,
     int *nfout2, float **arr2, float *tamp, int *nrow, int *ncol,
     float **ampo, int *ierr, char *name, char *pref);


void readhead_(int *sac, char *name, int *n, int *ierr);

void readdata_(int *sac, char *name, float *dt, float *delta,
     float *t0, float *seis, int *ierr);

void swapn(unsigned char *b, int N, int nn);


#endif /* !AFTAN_H */

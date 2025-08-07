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


void aftanipg_(double *PIover4, int *n, double *seis, double *t0, double *dt,
    double *delta, double *vmin, double *vmax, double *tmin, double *tmax,
    double *tresh, double *ffact, double *perc, int *npoints,
    double *taperl, int *nfin, double *snr, double *fmatch,
    int *npred, double **pred,
    int *nprpv, double *prpvper, double *prpvvel,
    int *nfout1, double **arr1, int *nfout2, double *arr2,
    double *tamp, int *nrow, int *ncol, double **ampo, int *ierr);

void aftanpg_(double *PIover4, int *n, double *seis, double *t0, double *dt,
    double *delta, double *vmin, double *vmax, double *tmin, double *tmax,
    double *tresh, double *ffact, double *perc, int *npoints,
    double *taperl, int *nfin, double *snr,
    int *nprpv, double *prpvper, double *prpvvel,
    int *nfout1, double **arr1, int *nfout2, double *arr2,
    double *tamp, int *nrow, int *ncol, double **ampo, int *ierr);

/*void printres_(double *dt, double *delta, int *nfout1, double **arr1,
 *     int *nfout2, double **arr2, double *tamp, int *nrow, int *ncol,
 *     double **ampo, int *ierr, char *name, char *pref);*/
void printres_(double *dt, double *delta, int *nfout1, double *arr1,
     int *nfout2, double *arr2, double *tamp, int *nrow, int *ncol,
     double *ampo, int *ierr, char *name, char *pref);


void readhead_(int *sac, char *name, int *n, int *ierr);

void readdata_(int *sac, char *name, double *dt, double *delta,
     double *t0, float *seis, int *ierr);

void swapn(unsigned char *b, int N, int nn);


#endif /* !AFTAN_H */

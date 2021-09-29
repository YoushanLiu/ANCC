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

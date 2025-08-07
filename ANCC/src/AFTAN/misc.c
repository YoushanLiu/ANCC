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
/*
 * Read SAC header and seismogram and make output
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include "aftan.h"
#include "cucss2sac.h"


void readhead_(int *sac, char *name, int *n, int *ierr)
{

    int i, nn, iswap, *nsam;
    int i1 = 0, i2 = 0, idelta = 0;

    FILE *fd;
    SAC_HD sach; /* the SAC header     */

    float dumm, t0, dt;

    float dum[2];

    char *p = 0, ch = 0;


    *ierr = 0;
    for (i = 0; i < 160; i++)
    {
        if (name[i] == '\0' || name[i] == ' ')
        {
            p = &name[i];
            ch = *p;
            break;
        }
    }
    *p = 0;

    iswap = 0;
    if ((fd = fopen(name, "r")) == NULL)
    {
        printf("Can not find file %s.\n", name);
        exit(1);
    }
    *p = ch;

    if (*sac)
    {
        /* The header */
        fread(&sach, sizeof(SAC_HD), 1, fd);
        if ((sach.internal4 > 100) || (sach.internal4 < -100))
            iswap = 1;
        if (iswap)
            swapn((unsigned char *)&sach, (int)(sizeof(float)), 158);
        dt = sach.delta;
        nsam = &sach.npts;
        t0 = sach.b;
        *n = *nsam;
        /* make symmetric cross-correlaton */
        nn = (1.0 - t0 / dt) + 0.5;
        if ((2 * nn - 1) != *n)
        {
            /* lag is not consistent with # of samples */
            printf("Wrong SAC header n=%d lag=%d\n", *n, (int)(-t0));
            *ierr = 1;
            return;
        }
        *n = nn;
    }
    fclose(fd);

}


void readdata_(int *sac, char *name, double *dt, double *delta,
               double *t0, float *seis, int *ierr)
{

    int i, n, nn, iswap, *nsam;
    int i1 = 0, i2 = 0, idelta = 0;

    FILE *fd;
    SAC_HD sach; /* the SAC header     */

    float dumm;
    float dum[2];
    float *tr;
    char *p = 0, ch = 0;


    *ierr = 0;
    for (i = 0; i < 160; i++)
    {
        if (name[i] == '\0' || name[i] == ' ')
        {
            p = &name[i];
            ch = *p;
            break;
        }
    }
    *p = 0;

    iswap = 0;
    if ((fd = fopen(name, "r")) == NULL)
    {
        printf("Can not find file %s.\n", name);
        exit(1);
    }
    *p = ch;

    if (*sac)
    {
        /* The header */
        fread(&sach, sizeof(SAC_HD), 1, fd);
        if ((sach.internal4 > 100) || (sach.internal4 < -100))
            iswap = 1;
        if (iswap)
            swapn((unsigned char *)&sach, (int)(sizeof(float)), 158);
        *dt = sach.delta;
        *delta = sach.dist;
        nsam = &sach.npts;
        *t0 = sach.b;
        /* The body */
        tr = (float *)malloc(*nsam * sizeof(float));
        fread(tr, sizeof(float), *nsam, fd);
        if (iswap)
            swapn((unsigned char *)tr, (int)(sizeof(float)), *nsam);
        n = *nsam;
        /* make symmetric cross-correlaton */
        nn = (1.0 - *t0 / *dt) + 0.5;
        if ((2 * nn - 1) != n)
        {
            /* lag is not consistent with # of samples */
            printf("Wrong SAC header n=%d lag=%d\n", n, (int)(-*t0));
            *ierr = 1;
            return;
        }
        for (i = 0; i < nn; i++)
            tr[n-i-1] = 0.50f*(tr[i] + tr[n-i-1]);
        for (i = 0; i < nn; i++)
            seis[i] = tr[nn-1+i];
        seis[0] *= 0.5f;
        //*n = nn;
        *t0 = 0.0;
        free(tr);
    }
    else
    {
        /*
         * Read test data from ascii file
        */
        i = 0;
        while ((nn = fscanf(fd, "%f %f", &dumm, &seis[i])) == 2)
        {
            if (i < 2)
                dum[i] = dumm;
            i++;
        }
        /* extract delta from file name */
        p = strrchr(name, '/');
        sscanf(p + 1, "proba_%d_%d_%d.SAC", &i1, &i2, &idelta);
        *delta = idelta;
        *dt = dum[1] - dum[0];
        //n = i;
        *t0 = 0.3;
    }
    fclose(fd);

}



/*
 * print completion result
 */
/*void printres_(double *xdt, double *delta, int *xnfout1, double **arr1,
 *               int *xnfout2, double **arr2, double *xtamp, int *xnrow,
 *               int *xncol, double **ampo, int *xierr, char *xname, char *xpref)*/
void printres_(double *xdt, double *delta, int *xnfout1, double *arr1,
               int *xnfout2, double *arr2, double *xtamp, int *xnrow, int *xncol, 
               double *ampo, int *xierr, char *xname, char *xpref)
{

    int i, j, k, ip, ierr;
    int nfout1, nfout2, nrow, ncol;

    FILE *out;

    float dt, tamp;


    char *p = 0, ch = 0, name2[160], name[160], pref[5];


    strncpy(pref, xpref, 2);
    pref[2] = '\0';
    for (i = 0; i < 160; i++)
    {
        if (xname[i] == '\0' || xname[i] == ' ')
        {
            p = &xname[i];
            ch = *p;
            break;
        }
    }

    *p = 0;
    strcpy(name, xname);
    *p = ch;
    dt = *xdt;
    nfout1 = *xnfout1;
    nfout2 = *xnfout2;
    tamp = *xtamp;
    nrow = *xnrow;
    ncol = *xncol;
    ierr = *xierr;

    /* write results to hard drive
       file ...DISP.0 contains preliminary result */

    strcpy(name2, name);
    strcat(name2, pref);
    strcat(name2, "_DISP.0");

    if ((out = fopen(name2, "w")) == NULL)
    {
        printf("Can not open file %s.\n", name2);
        exit(1);
    }
    for (i = 0; i < nfout1; i++)
    {
        //fprintf(out, "%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf\n",
        //    i, arr1[i][0], arr1[i][1], arr1[i][2], arr1[i][3],
        //    arr1[i][4], arr1[i][5], arr1[i][6]);
        k = i*8;
        fprintf(out, "%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf\n",
            i, arr1[k], arr1[k+1], arr1[k+2], arr1[k+3], arr1[k+4], arr1[k+5], arr1[k+6]);
    }
    fclose(out);


    /* file ...DISP.1 includes final results */
    if (nfout2 != 0)
    {
        strcpy(name2, name);
        strcat(name2, pref);
        strcat(name2, "_DISP.1");
        if ((out = fopen(name2, "w")) == NULL)
        {
            printf("Can not open file %s.\n", name2);
            exit(1);
        }
        for (i = 0; i < nfout2; i++)
        {
            //fprintf(out, "%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf\n",
            //    i, arr2[i][0], arr2[i][1], arr2[i][2], arr2[i][3],
            //    arr2[i][4], arr2[i][5]);
            k = i*7;
            fprintf(out, "%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf\n",
                i, arr2[k], arr2[k+1], arr2[k+2], arr2[k+3], arr2[k+4], arr2[k+5]);
        }
        fclose(out);
    }

    /* Output amplitude array into file on hard drive */
    strcpy(name2, name);
    strcat(name2, pref);
    strcat(name2, "_AMP");
    if ((out = fopen(name2, "w")) == NULL)
    {
        printf("Can not open file %s.\n", name2);
        exit(1);
    }

    fprintf(out, "%5d %5d %10.4lf %12.5lf\n", nrow, ncol, dt, *delta);
    for (i = 0; i < nrow; ++i)
    {
        ip =i*8 + 1;
        k = i*ncol;
        for (j = 0; j < ncol; ++j)
        	{
            	//fprintf(out, "%5d %8.3lf %15.6e\n", arr2[i][0], *delta / (tamp + j * dt), ampo[i][j]);
            	fprintf(out, "%10.4lf %8.3lf %15.6e\n", arr1[ip], *delta / (tamp + j * dt), ampo[k+j]);
        	}
    }
    fclose(out);

}

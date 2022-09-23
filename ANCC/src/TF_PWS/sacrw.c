/*
This file is part of ANCC.
TF_PWS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
TF_PWS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#define MAIN
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mysac.h"


#define PI 3.141592653589793
#define R 6371.0


void read_sac(char *fname, float *wf, SACHead *shd, long Nmax)
{

    FILE *fsac;

    fsac = fopen(fname, "rb");

    if (!fsac)
    {
        printf("The %s does not exit ! \n", fname);
        exit(0);
    }

    if (!shd) shd = &SAC_HEADER;
    fread(shd, sizeof(SACHead), 1, fsac);
    if (shd->npts > Nmax)
    {
        shd->npts = Nmax;
    }

    fread(wf, sizeof(float), (int)(shd->npts), fsac);

    fclose(fsac);

    //return shd;

}


void write_sac(char *fname, float *wf, SACHead *shd)
{

    int i, n;

    FILE *fsac;

    fsac = fopen(fname, "wb");

    if ((fsac = fopen(fname, "wb")) == NULL)
    {
        printf("Could not open sac file to write \n");
        exit(1);
    }

    if (!shd)
    {
        shd = &SAC_HEADER;
    }


    shd->iftype = (long)ITIME;
    shd->leven = (long)TRUE;

    shd->lovrok = (long)TRUE;
    shd->internal4 = 6L;




    shd->depmin = wf[0];
    shd->depmax = wf[0];

    for (i = 0; i < shd->npts; i++)
    {
        if (shd->depmin > wf[i]) shd->depmin = wf[i];
        if (shd->depmax < wf[i]) shd->depmax = wf[i];
    }

    // compute distance
    geodist(&SACHead);


    fwrite(shd, sizeof(SACHead), 1, fsac);

    n = fwrite(wf, sizeof(float), (int)(shd->npts), fsac);
    if (n != shd->npts)
    {
        printf("Error occur when writing the SAC ! ");
        printf("n:%d   shd -> npts:%d \n", n, (int)(shd->npts));
    }

    fclose(fsac);

}


void geodist(SACHead *shd)
{
    double stla, stlo, evla, evlo, deg2rad;
    //float PI = 3.141592653589793;

    deg2rad = PI / 180.0;

    stla = shd->stla * deg2rad;
    stlo = shd->stlo * deg2rad;
    evla = shd->evla * deg2rad;
    evlo = shd->evlo * deg2rad;

    shd->dist = (float)(R * acos(sin(stla)*sin(evla) + cos(stla)*cos(evla)*cos(stlo - evlo)));

}

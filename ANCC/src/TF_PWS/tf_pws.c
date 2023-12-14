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

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include "mysac.h"
//#include <omp.h>
#include <unistd.h>


#define E   1.e-4
#define EPS 1.e-16
#define MAINE
#define Nmax 131072
#define USAGE "[-B/-E: minimum or maximum period for group or phase velocity] \
               [-F: input file] [-O: output file] [-W: the exponential term of the weight]\n [-I: 1 or 2 inverse algorithm]"


#define  max(x,y)    ((x) > (y) ? (x) : (y))
#define  min(x,y)    ((x) < (y) ? (x) : (y))



// compile: gcc tf-PWS.c S-stform.c sacrw.c -lfftw3 -lm -o tf-PWS

/*
 * 	This code is written for tf-PWS based on Stockwell et al. 1996
 * 	USAGE:
 * 	    default paramenters:
 * 	        A.   Normalization and smooth of the weight funciton is turned off, if you want to turn on them,
 *               you can appendix -N1 and -S10 in the tail tf-PWS;
 *          B.   The default weight is set to 1, if you want use other value(for example n), you can appendix
 *               -Wn in the tail of tf-PWS;
 *          C.   The default inverse algorithm is defined by stockwell et al. 1996, and if you want to use the
 *               inverse algorithm proposed by Schimmel and Gallart 2005, you can appendix -I2 in the tail of
 *               tf-PWS.
 *          D.   The default frequent band is set from 0-Nyquist Hz, if you want to sepcifiy in some frequency
 *               band (for example from 10s to 100s), you can appendix -B10 -E100s in the tail.
 *          E.   This code support the parallel computing, and the default CPU number is set to 4, you can change
 *               that value too.
 *
 *        Examples:
 * 	        (1) ls * | ./tf-PWS
 * 	       	        to perform the line and tf-PWS of all the input files, and all other paramenter
 * 	                are set as default values;
 * 	        (2) ls * > input.dat
 * 	                ./tf-PWS -Finput.dat -B10  -E100  -N1 -S10 -I1 -W1
 *	                 In this case, the interested band is 10~100s but not 0-Nyquist. By this mean, you can reduce the
 *	                 calcaulation time, and filter the data in band 10~100s.
 *
 * 	Written by Guoliang Li, guoleonlee@gmail.com, Modified by Youshan Liu, at Institute of Geology and Geophysics,
 *           Chinese Academy of Sciences.
 * 	Finish data: 2016.01.10
 * 	 recitified 2016.08.09  Adding parallel computing part, and the number of parallel CPU is 4.
 * 	 recitified 2016.09.09  Fix the bug of storage overflowing.
 * 	 modified on 2016.09.20 adding the second inverse algorithm proposed by Schimmel and Gallart 2005;
 */


int getopt();
extern char *optarg;
extern int   optind;
void Strans(int len, int lo, int hi, double df, float *data, float *result);
void iStrans(int len, int lo, int hi, float *data, float *result);
void iStrans2(int len, int lo, int hi, double df, float *data, float *result);


int main(int arg, char *arv[])
{

    int i, j, n, m, itag, freq, nfft;
    int ispws, I, flag = 0, npow = 1;
    int smooth, normalize, user1;
    int lo, hi, lag, N, nstack, nq;
    int npts, nfreq, nst, ioffset;

    float MinP, MaxP, weight;
    float realpart, imagpart;
    float b, coeff, fmax, ampl;

    float *wf, *stack, *st;
    float **real, **imag;

    double df, dt, PI = 4.0*atan(1.0);
    double f1, f2, f3, f4;

    char infn[200], sacfn[200], outfn[200], outfile[200];

    // sac header
    SACHead shd;

    FILE *fp;



    if (arg < 1)
    {
        printf(USAGE);
        exit(0);
    }


    // Set the default parameters
    strcpy(outfile, "untitled");

    weight = 1.0;
    smooth = 0;
    normalize = 0;
    MinP = 0.0;
    MaxP = 0.0;
    f1 = 0;
    f4 = 0;
    I = 1;
    itag = 0;
    ispws = 0;


    while ((lag = getopt(arg, arv, "B:E:F:O:W:S:N:I:P:h")) != (-1))
    {
        switch (lag)
        {
        case 'B':
            MinP = atof(optarg);
            break;
        case 'E':
            MaxP = atof(optarg);
            break;
        case 'F':
            strcpy(infn, optarg);
            itag = 1;
            break;
        case 'O':
            strcpy(outfile, optarg);
            break;
        case 'W':
            weight = atof(optarg);
            break;
        case 'S':
            smooth = atoi(optarg);
            break;
        case 'N':
            normalize = atoi(optarg);
            break;
        case 'I':
            I = atoi(optarg);
            break;
        case 'P':
            ispws = atoi(optarg);
            break;
        case 'h':
            printf(USAGE);
            exit(0);
        }
    }


    if (MaxP < MinP)
    {
        printf("The input parameters are wrong set that should be MaxP < MinP, please check ! \n");
        exit(2);
    }

    if ((MaxP - MinP) > 1.e-4)
    {
        f1 = 1.0 / MaxP;
        f4 = 1.0 / MinP;
    }


    wf = (float *)malloc(Nmax * sizeof(float));


    // open the input file
    if (1 == itag)
    {
        if ((fp = fopen(infn, "r")) == NULL)
        {
            printf("The file %s cannot be read ! \n", infn);
            exit(1);
        }
    }
    else
    {
        fp = stdin;
    }



    // calculate the weights and the linear stack
    user1 = 0;
    nstack = 0;


    // only for linear stack
    if (ispws < 1)
    {

        for (i = 0; ; i++)
        {
            if (fscanf(fp, "%s", sacfn) == EOF)
            {
                break;
            }

            read_sac(sacfn, wf, &shd, Nmax);

            if (1 != flag)
            {
                b = shd.b;
                dt = shd.delta;
                npts = shd.npts;

                stack = (float *)malloc(npts * sizeof(float));
                memset(stack, 0, sizeof(float)*npts);

                flag = 1;
            }

            for (n = 0; n < npts; n++)
            {
                stack[n] = stack[n] + wf[n];
            }

            nstack = nstack + 1;
        }

        if (0 == nstack)
        {
			free(wf);
			free(stack);
            exit(0);
        }


        // normalize the linear stack result
        if (nstack > 1)
        {
            for (n = 0; n < npts; n++)
                stack[n] = stack[n] / (float)nstack;
        }


        shd.user0 = nstack;

        sprintf(outfn, "%s%s", outfile, "_ls.SAC");
        write_sac(outfn, stack, &shd);

		free(wf);
		free(stack);
        exit(0);
 
    }


    // only for phase-weighted stacks
    for (i = 0; ; i++)
    {

        if (EOF == fscanf(fp, "%s", sacfn))
        {
            break;
        }

        if (0 != access(sacfn, F_OK)) continue;

        read_sac(sacfn, wf, &shd, Nmax);


        // allocate memory and initialize parameters according to the npts of input series
        if (1 != flag)
        {

            b = shd.b;
            dt = shd.delta;
            npts = shd.npts;

            // calculate the points for fft
            nfft = pow(2, ceil(log((double)npts) / log(2.0L)));
            nq = (int)(nfft / 2);
            df = 1.0L / (nfft*dt);


            if (nfft > Nmax)
            {
                printf("The nfft=%d is greater than the predefined Nmax=%d \n", nfft, Nmax);
                exit(2);
            }


            // Set lo (low frequency) and hi (high frequency) for ST and IST
            lo = (int)(f1 / df);
            hi = min(ceil(f4 / df), nq);

            if ((0 == lo) && (0 == hi))
            {
                lo = 0;
                hi = nq;
                //hi = (int)(0.5 / (dt*df));  // 1.0/(2.0*dt) is the Nyquist frequency
                f4 = hi * df;
            }


            nfreq = hi - lo + 1;
            N = nfreq * nfft;


            // allocate linear stack array and initialize to zero
            stack = (float *)malloc(nfft * sizeof(float));
            memset(stack, 0, sizeof(float)*nfft);


            // allocate S-stform array and initialize to zero
            nst = 2 * N;
            st = (float *)calloc(nst, sizeof(float));
            memset(st, 0, sizeof(float)*nst);


            // allocate the weight arrays and initialize to zero
            real = (float **)malloc(nfreq * sizeof(float *));
            imag = (float **)malloc(nfreq * sizeof(float *));
            for (n = 0; n < nfreq; n++)
            {
                real[n] = (float *)malloc(nfft * sizeof(float));
                imag[n] = (float *)malloc(nfft * sizeof(float));
                memset(real[n], 0, sizeof(float)*nfft);
                memset(imag[n], 0, sizeof(float)*nfft);
            }

            flag = 1;

        }

        npts = (npts == shd.npts ? npts : 0);
        if (0 == npts)
        {
            exit(3);
        }

        // initialize the trailing to zero
        for (n = npts; n < min(nfft, Nmax); n++)
            wf[n] = 0.0;


        // linear stack
        for (n = 0; n < npts; n++)
        {
            stack[n] = stack[n] + wf[n];
        }


        // pre-filter
        //filter4_(&f1, &f2, &f3, &f4, &npow, &dt, &nfft, stack, stack);


        // Forward S-stransform
        Strans(nfft, lo, hi, df, &wf[0], &st[0]);


        //#pragma omp parallel for private(i, m, n, ioffset, ampl, coeff, realpart, imagpart) //num_threads(4)
        for (m = 0; m < nfreq; m++)
        {

            ioffset = m * 2*nfft;
            for (n = 0; n < nfft; n++)
            {

                i = 2 * n + ioffset;
                realpart = st[i    ];
                imagpart = st[i + 1];
                ampl = sqrt(realpart*realpart + imagpart*imagpart);
                coeff = 1.0 / max(ampl, EPS);
                //realpart = coeff * st[i    ];
                //imagpart = coeff * st[i + 1];
                realpart = coeff * realpart;
                imagpart = coeff * imagpart;

                // Prevent the emergence of the 0/0=Nan
                if ((fabs(realpart) > 1.0) || (fabs(imagpart) > 1.0))
                {
                    real[m][n] = real[m][n] + 0.0;
                    imag[m][n] = imag[m][n] + 0.0;
                }
                else
                {
                    real[m][n] = real[m][n] + realpart;
                    imag[m][n] = imag[m][n] + imagpart;
                }

            }

        }

        //// The variable user1 is used to store the stack times of tf-PWS, and finally are written in the SAC head as shd.user1
        //if (shd.user1 <= 0)
        //{
        //    user1 = user1 + 1;
        //}
        //else
        //{
        //    user1 = user1 + shd.user1;
        //}

        nstack++;

    }
    fclose(fp);



    if (0 == nstack)
    {
        exit(0);
    }


    // normalize the linear stack result
    if (nstack > 1)
    {
        for (n = 0; n < npts; n++)
            stack[n] = stack[n] / (float)nstack;
        for (n = npts; n < nfft; n++)
			stack[n] = 0.0;
    }


    // write linear stack result
    sprintf(outfn, "%s%s", outfile, "_ls.SAC");
    write_sac(outfn, stack, &shd);


    // calculate the spectrum of the linear stack result
    Strans(nfft, lo, hi, df, &stack[0], &st[0]);


    // calculate the weights
    //#pragma omp parallel for private(i, n, m, fmax, realpart, imagpart) //num_threads(4)
    for (i = 0; i < nfreq; i++)
    {

        fmax = 0.0;
        for (n = 0; n < nfft; n++)
        {

            // compute phase weight
            realpart = real[i][n] / (float)nstack;
            imagpart = imag[i][n] / (float)nstack;
            real[i][n] = pow(sqrt(realpart * realpart + imagpart * imagpart), weight);
            //real[i][n] = pow((realpart * realpart + imagpart * imagpart), 0.50*weight);

            if (real[i][n] > fmax)
            {
                fmax = real[i][n];
            }

        }

        // normalize along time axis
        if (1 == normalize)
        {
            for (n = 0; n < nfft; n++)
            {
                real[i][n] = real[i][n] / fmax;
            }
        }

    }

    // smooth along the time axis
    if (smooth > 1)
    {

        //smooth = 5;

        //#pragma omp parallel for private(i, n, m, ampl, lag) //num_threads(4)
        for (n = 0; n < nfreq; n++)
        {

            for (i = 0; i < nfft; i++)
            {

                lag = 0;
                ampl = 0.0;
                for (m = -smooth; m <= smooth; m++)
                {
                    if ((i + m >= 0) && (i + m < nfft))
                    {
                        ampl = ampl + real[n][i + m];
                        lag++;
                    }
                }

                if (lag > 0)
                {
                    imag[n][i] = ampl / (float)lag;
                }
                else
                {
                    imag[n][i] = real[n][i];
                }

                //if (fmax < imag[n][i])
                //{
                //    fmax = imag[n][i];
                //}

            }

        }

        for (n = 0; n < nfreq; n++)
        {
            for (i = 0; i < nfft; i++)
            {
                real[n][i] = imag[n][i];
            }
        }

    }


    //#pragma omp parallel for private(i, j, n, ioffset) //num_threads(4)
    for (i = 0; i < nfreq; i++)
    {

        ioffset = i * 2*nfft;
        for (n = 0; n < nfft; n++)
        {

            j = 2 * n + ioffset;

            // apply the weight to the spectrum of the linear stacking result
            st[j    ] = st[j    ] * real[i][n];
            st[j + 1] = st[j + 1] * real[i][n];

        }

    }

    shd.user0 = nstack;


    /* Perform the inverse tf-PWS: The inverse algorithm of iStrans is based on Stockwell er al. 1996,
       while that of iStrans2 is based on Schimmel et al. 2005, but we recommond the first inverse algorithm */

    // calculate inverse S-stransform
    if (I == 1)
        iStrans(nfft, lo, hi, &st[0], &stack[0]);

    if (I == 2)
        iStrans2(nfft, lo, hi, df, &st[0], &stack[0]);


    // filter4_(&f1, &f2, &f3, &f4, &npow, &dt, &nfft, stack, stack);

    sprintf(outfn, "%s%s", outfile, "_pws.SAC");
    write_sac(outfn, stack, &shd);


    // Free memory
    free(st);
    free(wf);
    free(stack);
    for (i = 0; i < nfreq; i++)
    {
        free(real[i]);
        free(imag[i]);
    }
    free(real);
    free(imag);


    return 0;

}

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
//#include <omp.h>

char *Wisfile = NULL;
char *Wistemplate = "%s/.fftwis";

#define WISLEN 8


void set_wisfile(void)
{

    char *home;

    if (Wisfile) return;

    home = getenv("HOME");

    Wisfile = (char *)malloc(strlen(home) + WISLEN + 1);
    sprintf(Wisfile, Wistemplate, home);

}

/* Convert frequencies in Hz into rows of the ST, given sampling rate and
   length. */
int st_freq(double f, int len, double srate)
{
    return floor(f * len / srate + .5);
}


static double gauss(int n, int m);


/* Stockwell transform of the real array data. The len argument is the
   number of time points, and it need not be a power of two. The lo and hi
   arguments specify the range of frequencies to return, in Hz. If they are
   both zero, they default to lo = 0 and hi = len / 2. The result is
   returned in the complex array result, which must be preallocated, with
   n rows and len columns, where n is hi - lo + 1. For the default values of
   lo and hi, n is len / 2 + 1. */
void Strans(int len, int lo, int hi, double df, float *data, float *result)
{

    int i, j, k, n, nq, l2;
    int nsqd, num, ioffset;

    float s, mean;
    float twopi, freq, scale;

    float *g;

    FILE *wisdom;

    static int planlen = 0;


    // fftwf_plan p1, p2;
    //fftwf_complex *cG;
    static fftwf_plan p1, p2;
    static fftwf_complex *h, *H, *G;


    nq = (int)(len / 2);
    if ((0 == lo) && (0 == hi)) {
        hi = nq;
    }


    /* Keep the arrays and plans around from last time, since this
        is a very common case. Reallocate them if they change. */
    if ((len != planlen) && (planlen > 0)) {
        fftwf_destroy_plan(p1);
        fftwf_destroy_plan(p2);
        fftwf_free(h);
        fftwf_free(H);
        fftwf_free(G);
        planlen = 0;
    }


    if (0 == planlen) {

        planlen = len;
        h = fftwf_malloc(sizeof(fftwf_complex) * len);
        H = fftwf_malloc(sizeof(fftwf_complex) * len);
        G = fftwf_malloc(sizeof(fftwf_complex) * len);

        /* Get any accumulated wisdom. */
        set_wisfile();

        wisdom = fopen(Wisfile, "r");
        if (wisdom) {
            fftwf_import_wisdom_from_file(wisdom);
            fclose(wisdom);
        }

        /* Set up the fftw plans. */
        p1 = fftwf_plan_dft_1d(len, h, H, FFTW_FORWARD, FFTW_MEASURE);
        p2 = fftwf_plan_dft_1d(len, G, h, FFTW_BACKWARD, FFTW_MEASURE);

        /* Save the wisdom. */
        wisdom = fopen(Wisfile, "w");
        if (wisdom) {
            fftwf_export_wisdom_to_file(wisdom);
            fclose(wisdom);
        }

    }


    /* Convert the input to complex. Also compute the mean. */
    s = 0.0;
    memset(h, 0, sizeof(fftwf_complex) * len);
    for (i = 0; i < len; i++) {
        h[i][0] = data[i];
        s += data[i];
    }
    s /= len;
    mean = s;

    /* FFT. */
    //fftwf_execute(p1); /* h -> H */
    fftwf_execute_dft(p1, h, H); /* h -> H */


    /* Hilbert transform. The upper half-circle gets multiplied by
       two, and the lower half-circle gets set to zero.
       The real axis is left alone. */
    l2 = (len + 1) / 2;
    for (i = 1; i < l2; i++) {
        H[i][0] *= 2.;
        H[i][1] *= 2.;
    }
    l2 = nq + 1;
    for (i = l2; i < len; i++) {
        H[i][0] = 0.;
        H[i][1] = 0.;
    }


    /* Fill in rows of the result. */

    /* The row for lo == 0 contains the mean. */
    /* Subsequent rows contain the inverse FFT of the spectrum
       multiplied with the FFT of scaled gaussians. */

    twopi = 2.0*M_PI*M_PI;
    //#pragma omp parallel for private(i, j, k, ioffset, n, nsqd, s, cG, g) //num_threads(4)
    #pragma omp parallel for private(i, j, k, ioffset, n, nsqd, s, G, g) //num_threads(4)
    for (n = lo; n <= hi; n++)
    {

        ioffset = (n - lo) * 2*len;

        nsqd = n * n;
        /* Scale the FFT of the gaussian. Negative frequencies
           wrap around. */

        /* The row for lo == 0 contains the mean. */
        /* Subsequent rows contain the inverse FFT of the spectrum
           multiplied with the FFT of scaled gaussians. */
        freq = n * df;

        /*
        if( freq > 0.025 ) scale=4.0;
        if( freq <=0.025 && freq>=0.02) scale=4.0-0.3*(1./freq-40);
        if( freq < 0.02 && freq>=0.01 ) scale=1.0-(1.0/freq-50)*0.01;
        if ( freq <0.01 ) scale=0.5;
        */

        /*Parameter scale controls the time-frequency resolution. Bigger value can help to
          increase the frequency resolution but decrease the time resolution, and vice versa
        */

        scale = 1.0;
        g = (float *)malloc(sizeof(float)*len);
        //cG = (fftwf_complex*)malloc(sizeof(fftwf_complex)*len);
        memset(G, 0, sizeof(fftwf_complex) * len);
        if (n == 0)
        {

            for (i = 0; i < len; i++)
            {
                j = 2 * i + ioffset;
                result[j    ] = mean;
                result[j + 1] = 0.0;
            }

        }
        else
        {

            g[0] = 1.0;

            //l2 = nq + 1;
            //for (i = 1; i < l2; i++)
            //{
            //    g[i] = g[len - i] = exp(-twopi * i * i * scale / nsqd);
            //}
            for (i = 1; i <= nq; i++)
            {
                g[i] = exp(-twopi * i * i * scale / nsqd);
            }
            for (i = nq + 1; i < len; i++)
            {
                g[i] = g[len - i];
            }

            for (i = 0; i < len; i++)
            {
                s = g[i];
                k = n + i;
                if (k >= len) k -= len;
                //cG[i][0] = s * H[k][0];
                //cG[i][1] = s * H[k][1];
                G[i][0] = s * H[k][0];
                G[i][1] = s * H[k][1];
            }

            /* Inverse FFT the result to get the next row. */
            //#pragma omp critical (setction1)
            {

				/*
                for (i = 0; i < len; i++)
                {
                    G[i][0] = cG[i][0];
                    G[i][1] = cG[i][1];
                }*/
                //fftwf_execute(p2); /* G -> h */
                fftwf_execute_dft(p2, G, h); /* G -> h */

                for (i = 0; i < len; i++)
                {

                    j = 2 * i + ioffset;
                    *(result + (j    )) = h[i][0] / (float)len;
                    *(result + (j + 1)) = h[i][1] / (float)len;

                }

            }
            /* Go to the next row. */

        }

        free(g);
        //fftwf_free(cG);

    }

    if (0 != planlen) {
    	fftwf_free(h);
    	fftwf_free(G);
    	fftwf_free(G);
	}
}

/* This is the Fourier Transform of a Gaussian. */
static double gauss(int n, int m)
{
    return exp(-2.0 * M_PI * M_PI * m * m / (n * n));
}



/* Inverse Stockwell transform, this inverse algorithm is based on Stockwell et al. 1996 */
void iStrans(int len, int lo, int hi, float *data, float *result)
{

    int i, j, n, nq;
    int ioffset, l2;

    float *p;

    FILE *wisdom;

    static int planlen = 0;
    static fftwf_plan p2;
    static fftwf_complex *h, *H;


#ifdef debuge
    FILE *fp, *fp1;

    if ((fp = fopen("new_SP.dat", "w")) == NULL)
    {
        printf("new_SP.dat cannot be writen!\n");
        exit(0);
    }
    if ((fp1 = fopen("signal.dat", "w")) == NULL)
    {
        printf("signal.dat cannot be writen!\n");
        exit(0);
    }
#endif


    nq = (int)(len / 2);
    /* Check for frequency defaults. */
    if ((lo == 0) && (hi == 0)) {
        hi = nq;
    }


    /* Keep the arrays and plans around from last time, since this
       is a very common case. Reallocate them if they change. */
    if ((len != planlen) && (planlen > 0)) {
        fftwf_destroy_plan(p2);
        fftwf_free(h);
        fftwf_free(H);
        planlen = 0;
    }


    if (0 == planlen) {

        planlen = len;
        h = fftwf_malloc(sizeof(fftwf_complex) * len);
        H = fftwf_malloc(sizeof(fftwf_complex) * len);

        /* Get any accumulated wisdom. */
        set_wisfile();
        wisdom = fopen(Wisfile, "r");
        if (wisdom) {
            fftwf_import_wisdom_from_file(wisdom);
            fclose(wisdom);
        }


        /* Set up the fftw plans. */
        p2 = fftwf_plan_dft_1d(len, H, h, FFTW_BACKWARD, FFTW_MEASURE);


        /* Save the wisdom. */
        wisdom = fopen(Wisfile, "w");
        if (wisdom) {
            fftwf_export_wisdom_to_file(wisdom);
            fclose(wisdom);
        }

    }


    /* Sum the complex array across time. */
    memset(H, 0, sizeof(fftwf_complex) * len);

    p = data;
    for (n = lo; n <= hi; n++)
    {

        for (i = 0; i < len; i++)
        {

#ifdef debuge
            fprintf(fp, "%8.5f ", *p);
#endif

            H[n][0] += *p++;

#ifdef debuge
            fprintf(fp, "%8.5f ", *p);
#endif

            H[n][1] += *p++;

        }

#ifdef debuge
        fprintf(fp, "\n");
#endif

    }

    /* Invert the Hilbert transform. */
    /*
    l2 = (len + 1) / 2;
    for (i = 1; i < l2; i++) {
        H[i][0] /= 2.;
        H[i][1] /= 2.;
    }
    l2 = len / 2 + 1;
    for (i = l2; i < len; i++) {
        H[i][0] = H[len - i][0];
        H[i][1] = -H[len - i][1];
    }
    for(i=0;i<hi;i++)
        printf("%d %7.5lf  %7.5lf %7.5lf  %7.5lf\n",i, H[i][0],H[i][1] ,H[len-1-i][0],H[len-1-i][1]);
    getchar();
    */
    /* Inverse FFT. */


    //fftwf_execute(p2); /* H -> h */
    fftwf_execute(p2, H, h); /* H -> h */
    fftwf_destroy_plan(p2);

    p = result;
    for (i = 0; i < len; i++) {

        *p++ = h[i][0] / (float)len;

#ifdef debuge
        fprintf(fp1, "%7.5lf ", h[i][0] / len);
#endif

    }

#ifdef debuge
    fclose(fp);
    fclose(fp1);
#endif

}


/* This does just the Hilbert transform. */
void hilbert(int len, float *data, float *result)
{
    int i, l2, nq;

    float *p;

    FILE *wisdom;

    static int planlen = 0;
    static fftwf_plan p1, p2;
    static fftwf_complex *h, *H;


    nq = (int)(len / 2);
    /* Keep the arrays and plans around from last time, since this
       is a very common case. Reallocate them if they change. */
    if ((len != planlen) && (planlen > 0)) {
        fftwf_destroy_plan(p1);
        fftwf_destroy_plan(p2);
        fftwf_free(h);
        fftwf_free(H);
        planlen = 0;
    }


    if (0 == planlen) {

        planlen = len;

        h = fftwf_malloc(sizeof(fftwf_complex) * len);
        H = fftwf_malloc(sizeof(fftwf_complex) * len);

        /* Get any accumulated wisdom. */
        set_wisfile();
        wisdom = fopen(Wisfile, "r");
        if (wisdom) {
            fftwf_import_wisdom_from_file(wisdom);
            fclose(wisdom);
        }

        /* Set up the fftw plans. */
        p1 = fftwf_plan_dft_1d(len, h, H, FFTW_FORWARD, FFTW_MEASURE);
        p2 = fftwf_plan_dft_1d(len, H, h, FFTW_BACKWARD, FFTW_MEASURE);

        /* Save the wisdom. */
        wisdom = fopen(Wisfile, "w");
        if (wisdom) {
            fftwf_export_wisdom_to_file(wisdom);
            fclose(wisdom);
        }

    }


    /* Convert the input to complex. */
    memset(h, 0, sizeof(fftwf_complex) * len);
    for (i = 0; i < len; i++) {
        h[i][0] = data[i];
    }


    /* FFT. */
    //fftwf_execute(p1); /* h -> H */
    fftwf_execute_dft(p1, h, H); /* h -> H */


    /* Hilbert transform. The upper half-circle gets multiplied by
       two, and the lower half-circle gets set to zero.
       The real axisis left alone. */

    l2 = (len + 1) / 2;
    for (i = 1; i < l2; i++) {
        H[i][0] *= 2.;
        H[i][1] *= 2.;
    }
    l2 = nq + 1;
    for (i = l2; i < len; i++) {
        H[i][0] = 0.;
        H[i][1] = 0.;
    }


    /* Inverse FFT. */
    //fftwf_execute(p2); /* H -> h */
    fftwf_execute_dft(p2, H, h); /* H -> h */


    /* Fill in the rows of the result. */
    p = result;
    for (i = 0; i < len; i++) {
        *p++ = h[i][0] / (float)len;
        *p++ = h[i][1] / (float)len;
    }

}


/* Inverse Stockwell transform, this inverse algorithm is based on Schimmel et al. 2005,
   and this inverse algorihm is only used for testing */
void iStrans2(int len, int lo, int hi, double df, float *data, float *result)
{

    int i, j, n, nq, l2;

    float *p;

    float out, mean;
    float coeff, scale;

    FILE *wisdom;

    static int planlen = 0;
    static fftwf_plan p2;
    static fftwf_complex *h, *H;


    /*
    FILE *fp;
    if((fp=fopen("new_SP.dat","w"))==NULL)
    {
            printf("new_SP.dat cannot be writen!\n");
            exit(0);
    }
    */

#ifdef debuge
    FILE *fp, *fp1;

    if ((fp = fopen("new_SP.dat", "w")) == NULL)
    {
        printf("new_SP.dat cannot be writen!\n");
        exit(0);
    }
    if ((fp1 = fopen("signal.dat", "w")) == NULL)
    {
        printf("signal.dat cannot be writen!\n");
        exit(0);
    }
#endif


    nq = (int)(len / 2);
    /* Check for frequency defaults. */
    if ((0 == lo) && (0 == hi)) {
        hi = nq;
    }


    /* Keep the arrays and plans around from last time, since this
       is a very common case. Reallocate them if they change. */
    if ((len != planlen) && (planlen > 0)) {
        fftwf_destroy_plan(p2);
        fftwf_free(h);
        fftwf_free(H);
        planlen = 0;
    }


    if (0 == planlen) {

        planlen = len;

        h = fftwf_malloc(sizeof(fftwf_complex) * len);
        H = fftwf_malloc(sizeof(fftwf_complex) * len);

        /* Get any accumulated wisdom. */
        set_wisfile();
        wisdom = fopen(Wisfile, "r");
        if (wisdom) {
            fftwf_import_wisdom_from_file(wisdom);
            fclose(wisdom);
        }

        /* Set up the fftw plans. */
        p2 = fftwf_plan_dft_1d(len, H, h, FFTW_BACKWARD, FFTW_MEASURE);

        /* Save the wisdom. */
        wisdom = fopen(Wisfile, "w");
        if (wisdom) {
            fftwf_export_wisdom_to_file(wisdom);
            fclose(wisdom);
        }

    }



    /* Sum the complex array across time. */
    scale = sqrt(2.0 * M_PI);
    memset(H, 0, sizeof(fftwf_complex) * len);


    p = data;
    for (i = 0; i < len; i++)
    {

        for (n = lo; n <= hi; n++)
        {

            j = 2*(i + (n - lo)*len);

            H[n][0] = *(p + j);
            H[n][1] = *(p + j + 1);
            if (0 == n)
            {
                mean = H[n][0];
                H[n][0] = 0;
                H[n][1] = 0;
            }
            if (n > 0)
            {
            	coeff = scale / ((float)n*df);
                H[n][0] = coeff * H[n][0];
                H[n][1] = coeff * H[n][1];
            }

        }

        //fftwf_execute(p2); /* H -> h */
        fftwf_execute(p2, H, h); /* H -> h */

        out = 0.0;
        for (l2 = 0; l2 < len; l2++)
        {
            out = out + h[l2][0];
        }

        /*
        if( i%40==10 )
        for(n=0;n<len;n++)
        {
            fprintf(fp,"%d  %d  %f\n",i,n,h[n][0]);
        }
        */

        *result++ = h[i][0] / (float)len + mean;

    }


    /* Invert the Hilbert transform. */
    /*
    l2 = (len + 1) / 2;
    for (i = 1; i < l2; i++) {
        H[i][0] /= 2.;
        H[i][1] /= 2.;
    }
    l2 = len / 2 + 1;
    for (i = l2; i < len; i++) {
        H[i][0] = H[len - i][0];
        H[i][1] = -H[len - i][1];
    }
    for(i=0;i<hi;i++)
        printf("%d %7.5lf  %7.5lf %7.5lf  %7.5lf\n",i, H[i][0],H[i][1] ,H[len-1-i][0],H[len-1-i][1]);
    getchar();
    */
    /* Inverse FFT. */


    fftwf_destroy_plan(p2);


#ifdef debuge
    fclose(fp);
    fclose(fp1);
#endif

}


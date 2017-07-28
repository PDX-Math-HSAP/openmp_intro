#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mm_malloc.h>

//#define OPTREPORT

#define NTHREADS 4

#ifndef OPTREPORT
double implementation0(double * a, double * b, int N);
double implementation1(double * a, double * b, int N);
double implementation2(double * a, double * b, int N);
double implementation3(double * a, double * b, int N);
double implementation4(double * a, double * b, int N);
double implementation5(double * a, double * b, int N);
double implementation6(double * a, double * b, int N);
double implementation7(double * a, double * b, int N);
double implementation8(double * a, double * b, int N);
#endif

int main(int argc, char *argv[])
{
    int N = 300000000;
    int nruns = 1;
    int max_threads;
    //int run, i;
    double average_time;
    double standard_deviation;
    double delta;
    double sec;
    double * a = (double*)_mm_malloc(N * sizeof(double), 16);
    double * b = (double*)_mm_malloc(N * sizeof(double), 16);
    //double * a = (double*)malloc(N * sizeof(double));
    //double * b = (double*)malloc(N * sizeof(double));

    //srand(time(NULL));

#ifdef _OPENMP
    max_threads = omp_get_max_threads();
    printf ("Warning: we can use %d threads only \n", max_threads);
    //int num_threads = omp_get_num_threads();
    //printf ("Warning: we can use %d threads only \n", num_threads);
#endif


#ifdef OPTREPORT
    for (int i = 0; i < N; ++i)
        a[i] = cos(i * 1.0 / 200);//rand() * 1.0 / RAND_MAX;
    for (int i = 0; i < N; ++i)
	b[i] = cos(i * 1.0 / 200);//rand() * 1.0 / RAND_MAX;

    // This doesn't work
    clock_t start = clock(), diff;

    double sum = 0;

    for (int i = 0; i < N; ++i)
	sum += a[i]*b[i];

    diff = clock() - start;

    sec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    // But this is getting vectorized? No, unlucky :(

    sum = 0;
    double tmp[4];

    start = clock();

    for (int i = 0; i < N / 4; ++i)
    {
	tmp[0] = a[4*N + 0] * b[4*N + 0];
	tmp[1] = a[4*N + 1] * b[4*N + 1];
	tmp[2] = a[4*N + 2] * b[4*N + 2];
	tmp[3] = a[4*N + 3] * b[4*N + 3];

	sum = tmp[0] + tmp[1] + tmp[2] + tmp[3];
    }

    diff = clock() - start;

    sec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

#else
    // running the loop, impl 0
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation0(a,b,N);
        printf ( "sec = %e \n", sec);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 0 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);

/*
    // running the loop, impl 1
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation1(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 1 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);

    // running the loop, impl 2
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation2(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 2 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);

    // running the loop, impl 3
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation3(a,b,N);
        printf ( "sec = %e \n", sec);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 3 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);
*/

/*
    // running the loop, impl 4
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation4(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 4 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);

    // running the loop, impl 5
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation5(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 5 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);

    // running the loop, impl 6
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation6(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 6 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);

    // running the loop, impl 7
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation7(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 7 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);
*/
    // running the loop, impl 8
    average_time = 0;
    standard_deviation = 0;
    for (int run = 0; run < nruns; ++run)
    {
        for (int i = 0; i < N; ++i)
	    a[i] = cos((i + 2.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        for (int i = 0; i < N; ++i)
	    b[i] = cos((i + 3.0*run)*1.0 / 200);//rand() * 1.0 / RAND_MAX;
        sec = implementation8(a,b,N);
        delta = sec - average_time;
        average_time += delta / (run + 1);
        standard_deviation += delta * (sec - average_time); 
    }
    standard_deviation /= nruns - 1;
    printf ("average_time for implement. 8 = %f milliseconds, disp. = %f \n", average_time, standard_deviation);
#endif

    _mm_free(a);
    _mm_free(b);
    //free(a);
    //free(b);
    return 0;
}

#ifndef OPTREPORT
double implementation0(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

    for (int i = 0; i < N; ++i)
	sum += a[i]*b[i];

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation1(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

#pragma omp parallel for num_threads(NTHREADS)
    for (int i = 0; i < N; ++i)
	sum += a[i]*b[i];

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}


double implementation2(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

#pragma omp parallel for num_threads(NTHREADS)
    for (int i = 0; i < N; ++i)
    {
#pragma omp atomic        
	sum += a[i]*b[i];
    }

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation3(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    int max_threads = omp_get_max_threads();

    double sum = 0;

#pragma omp parallel num_threads(NTHREADS) //max_threads)
{
    //printf ( "I am %d \n", omp_get_thread_num());

#pragma omp for reduction (+:sum) 
    for (int i = 0; i < N; ++i)
    {
	sum += a[i]*b[i];
    }
}
    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation4(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

#pragma omp parallel for reduction (+:sum) num_threads(8)
    for (int i = 0; i < N; ++i)
    {
	sum += a[i]*b[i];
    }

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation5(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

#pragma omp parallel for reduction (+:sum) shared(a,b,N) default(none) num_threads(NTHREADS)
    for (int i = 0; i < N; ++i)
    {
	sum += a[i]*b[i];
    }

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation6(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

#pragma omp parallel for reduction (+:sum) shared(a,b,N) default(none) schedule(dynamic) num_threads(NTHREADS)
    for (int i = 0; i < N; ++i)
    {
	sum += a[i]*b[i];
    }

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation7(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum = 0;

    int M = 32;
#pragma omp parallel for reduction(+:sum) shared(a,b,N) num_threads(NTHREADS)
    for (int i=0; i<N/M; i++)
    {
        for (int j=i*M; j<(i+1)*M; j++)
        {
            sum+=a[j]*b[j];
        }
    }
    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum);

    return sec;
}

double implementation8(double * a, double * b, int N)
{
    clock_t start = clock(), diff;

    double sum1 = 0, sum2 = 0;

#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
        {
            for (int i=0; i<N/2; i++)
            {
                 sum1+=a[i]*b[i];
            }
        }
#pragma omp section
        {
            for (int i=N/2; i<N; i++)
            {
                 sum2+=a[i]*b[i];
            }
        }
    }

    diff = clock() - start;

    double sec = (diff * 1.0) / CLOCKS_PER_SEC;
    //printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("sum = %e \n", sum1 + sum2);

    return sec;
}

#endif


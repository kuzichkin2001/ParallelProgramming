#include <stdio.h>
#include <omp.h>
#include <locale.h>
#include "Simps.h"
#include <math.h>
#include <time.h>
#ifdef _WIN32
#include <Windows.h>

void sleep(int seconds)
{
    Sleep(1000 * seconds);
}

#endif
#define _USE_MATH_DEFINES


void labwork_4_task_3();
void labwork_4_task_5();
void labwork_4_task_7();

const int NUMBER_OF_THREADS = 4;
const int N = (int)1e6;
omp_lock_t lock;

int main()
{
    // labwork_4_task_3();
    // labwork_4_task_5();
    labwork_4_task_7();
}

void labwork_4_task_3()
{
    int n;
#pragma omp parallel num_threads(NUMBER_OF_THREADS)
    {
#pragma omp critical
        {
            n = omp_get_thread_num();
            printf("Thread %d\n", n);
        }
    }
}

void labwork_4_task_5()
{
    int n;
    omp_init_lock(&lock);

#pragma omp parallel num_threads(NUMBER_OF_THREADS) private(n)
    {
        n = omp_get_thread_num();
        omp_set_lock(&lock);
        printf("Beginning of closed section, thread %d\n", n);
        sleep(1 + 1. / 2 * ((n + 1) % 2));
        printf("End of closed section, thread %d\n", n);
        omp_unset_lock(&lock);
    }
    omp_destroy_lock(&lock);
}

double F(double x)
{
    return 1.0 / (1.0 + x * x);
}

void labwork_4_task_7()
{
    double Tms = clock();
    double Intgr = Simps(0, 1000000, 1000000000, F);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    printf("Time=%f sec\n", Tms);
    printf("Intgr=%f\n", Intgr);
}
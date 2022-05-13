// Parallel Programming. Third labwork.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdio.h>
#include <omp.h>
#include <locale.h>

#ifdef _WIN32
#include <Windows.h>
void sleep(int secs)
{
    Sleep(1000 * secs);
}
#endif

const int NUMBER_OF_THREADS = 8;
void labwork_3_task_3();
void labwork_3_task_5();
void labwork_3_task_7();

int main(int argc, char* argv[])
{
    //labwork_3_task_3();
    //labwork_3_task_5();
    labwork_3_task_7();
}

void labwork_3_task_3()
{
    int currentNumOfThreads, currentThreadIdx;
#pragma omp parallel num_threads(NUMBER_OF_THREADS) 
    {
        currentNumOfThreads = omp_get_num_threads();
        currentThreadIdx = omp_get_thread_num();

        if (currentThreadIdx == 0) printf("Total number of threads: %d\n", currentNumOfThreads);
        else printf("Current thread index: %d\n", currentThreadIdx);
    }
}

void labwork_3_task_5()
{
    int i;
#pragma omp parallel num_threads(NUMBER_OF_THREADS) private(i)
    {
#pragma omp for schedule(static, 5)
        for (i = 0; i < NUMBER_OF_THREADS; ++i)
        {
            printf("Thread %d made iteration %d\n", omp_get_thread_num(), i);
            sleep(1);
        }
    }
}

void labwork_3_task_7()
{
    int number = 0;
#pragma omp parallel num_threads(NUMBER_OF_THREADS)
    {
#pragma omp sections
        {
#pragma omp section
            {
                number = 1;
                printf("smth, n: %d\n", number);
            }
#pragma omp section
            {
                number = 2;
                printf("smth2, n: %d\n", number);
            }
#pragma omp section
            {
                number = 3;
                printf("smth3, n: %d\n", number);
            }
        }
        printf("N in thread %d: %d\n", omp_get_thread_num(), number);
    }
    printf("Result n: %d\n", number);
}

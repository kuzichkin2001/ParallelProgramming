// Parallel programming, Second labwork.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdio.h>
#include <omp.h>
#include <locale.h>

void labwork_2_task_3();
void labwork_2_task_5();
void labwork_2_task_7();
const int NUMBER_OF_THREADS = 4;
const int N = 20;

int main(int argc, char* argv[])
{
    //labwork_2_task_3();
    //labwork_2_task_5();
    labwork_2_task_7();
}

void labwork_2_task_3()
{
#pragma omp parallel num_threads(NUMBER_OF_THREADS)
    {
        printf("Parallel section 1\n");
#pragma omp single nowait
        {
            printf("Single thread\n");
        }
        printf("Parallel section 2\n");
    }
}

void labwork_2_task_5()
{
    int n;
#pragma omp parallel num_threads(NUMBER_OF_THREADS) private(n)
    {
        n = 1;
#pragma omp master
        {
            n = 2;
        }
        printf("First n: %d\n", n);
#pragma omp barrier
#pragma omp master
        {
            n = 3;
        }
        printf("Second n: %d\n", n);
    }
}

void labwork_2_task_7()
{
    int arr[N];
    printf("Initialized array:\n");

    for (int i = 0; i < N; ++i)
    {
        arr[i] = 0;
        printf("%d ", arr[i]);
    }

    printf("\n");

#pragma omp parallel num_threads(N) shared(arr)
    {
        int idx = omp_get_thread_num();
        arr[idx] = idx;
    }

    printf("Correctly filled array:\n");
    for (int i = 0; i < N; ++i) printf("%d ", arr[i]);
}
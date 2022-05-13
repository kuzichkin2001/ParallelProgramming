// Parallel programming.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdio.h>
#include <omp.h>

void labwork_1_task_3();
void labwork_1_task_5();
void labwork_1_task_7();
const int NUMBER_OF_THREADS = 7;

int main(int argc, char* argv[])
{
    //labwork_1_task_3();
    //labwork_1_task_5();
    labwork_1_task_7();

    return 0;
}


// Labwork 1, Task 3
void labwork_1_task_3()
{
    printf("Section 1\n");

    omp_set_num_threads(NUMBER_OF_THREADS);
#pragma omp parallel
    {
        printf("Section 2\n");
    }
}

// Labwork 1, Task 5
void labwork_1_task_5() {
    omp_set_dynamic(0);

    omp_set_num_threads(NUMBER_OF_THREADS);

    #pragma omp parallel num_threads(3)
    {
        printf("Parallel section 1\n");
    }

    #pragma omp parallel
    {
        printf("Parallel section 2\n");
    }
}

// Labwork 1, Task 7
void labwork_1_task_7()
{
    printf("%d\n", omp_get_dynamic());

    omp_set_dynamic(1);

    printf("Dynamic determined number of threads: %d\n", omp_get_dynamic());

    #pragma omp parallel num_threads(3)
    {
        #pragma omp master
        {
            printf("Number of threads in parallel section: %d\n", omp_get_num_threads());
        }
    }
}

#include <iostream>
#include <cmath>
#include <time.h>
#include <omp.h>

#define N 1000

using namespace std;

double func(const double&, const double&, const double&);

int main(int* argc, char** argv)
{

    double h = 1.0 / N;
    double A = 0, B = 0, Tms = clock();

#pragma omp parallel for reduction(+: A)
    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            for (int k = 0; k < N - 1; ++k)
            {
                A += func(h * (i + 1. / 2.), h * (j + 1. / 2.), h * (k + 1. / 2.));
            }
        }
    }

    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;

    Tms = clock();
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (static)
    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            for (int k = 0; k < N - 1; ++k)
            {
                A += func(h * (i + 1. / 2.), h * (j + 1. / 2.), h * (k + 1. / 2.));
            }
        }
    }

    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time static = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;

    Tms = clock();
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (dynamic)
    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            for (int k = 0; k < N - 1; ++k)
            {
                A += func(h * (i + 1. / 2.), h * (j + 1. / 2.), h * (k + 1. / 2.));
            }
        }
    }

    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time dynamic = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;

    Tms = clock();
    A = 0;

#pragma omp parallel for reduction(+: B)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N - 1; j++) {
            for (int k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
                B += func(h * i, h * j, h * k);
        }
    }
    B *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time = " << Tms << " sec" << endl;
    cout << "B = " << B << endl;
    B = 0;
    Tms = clock();

#pragma omp parallel for reduction(+: B) schedule (static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N - 1; j++) {
            for (int k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
                B += func(h * i, h * j, h * k);
        }
    }
    B *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time static = " << Tms << " sec" << endl;
    cout << "B = " << B << endl;
    B = 0;
    Tms = clock();

#pragma omp parallel for reduction(+: B) schedule (dynamic)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N - 1; j++) {
            for (int k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
                B += func(h * i, h * j, h * k);
        }
    }
    B *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time dynamic = " << Tms << " sec" << endl;
    cout << "B = " << B << endl;
    B = 0;
    Tms = clock();

    return 0;
}

double func(const double& x, const double& y, const double& z)
{
    return cos(x - pow(y, 2)) / (2 + exp(-pow(z, 2)) * sin(pow(x, 2) + pow(y, 2) - z));
}

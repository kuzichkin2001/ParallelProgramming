#include "Simps.h"
#include <stdio.h>
#include <omp.h>
#include <locale.h>


double Simps(double a, double b, int N, double Func(double))
{
	double h = (b - a) / (2 * N);
	int k;
	double S1 = 0, S2 = 0, Tmp;

#pragma omp parallel num_threads(4) shared(a, h) private(k, Tmp) \
	reduction(+:S1) reduction(+:S2)
	{
		for (k = 1; k < N; ++k)
		{
			Tmp = a + (2 * k - 1) * h;
			S1 += Func(Tmp);
			S2 += Func(Tmp + h);
		}

		S1 += Func(b - h);
	}

	return h * (Func(a)
		+ Func(b)
		+ 4.0 * S1
		+ 2.0 * S2) / 3.0;
}
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <cmath>
#include "MPI_Simps.h"

#define NN 1000
#define NNN 3

// labwork #6
void labwork_6_task_3();
void labwork_6_task_5();
void labwork_6_task_6();

// labwork #7
void labwork_7_task_3();
void labwork_7_task_5();
void labwork_7_task_6();

// labwork #8
void labwork_8_task_3();
void labwork_8_task_5();
void labwork_8_task_7();

int main(int argc, char* argv[])
{
    //labwork_6_task_3();
    //labwork_6_task_5();
    //labwork_6_task_6();

    //labwork_7_task_3();
    //labwork_7_task_5();
    //labwork_7_task_6();

    //labwork_8_task_3();
    labwork_8_task_5();
    //labwork_8_task_7();

    return 0;
}

void labwork_6_task_3()
{
    MPI_Init(NULL, NULL);
    int Rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    if (Rank == 0)
    {
        double Val = 7.77;
        int BufSize = sizeof(double) + MPI_BSEND_OVERHEAD;
        void* Buff = malloc(BufSize);

        MPI_Buffer_attach(Buff, BufSize);
        MPI_Bsend(&Val, 1, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD);
        MPI_Buffer_detach(Buff, &BufSize);

        if (!Buff) free(Buff);
        std::cout << "Process 0 has sent the data" << std::endl;
    }
    else
    {
        if (Rank == 1)
        {
            double Val = 0;
            MPI_Status Status;
            MPI_Recv(&Val, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &Status);
            std::cout << "Process 1 has received the data from process 0. Val = " << Val << std::endl;
        }
        else std::cout << "Process " << Rank << " does not participate in message exchanging." << std::endl;
    }

    MPI_Finalize();
}

void labwork_6_task_5()
{
    MPI_Init(NULL, NULL);

    int Size;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);

    if (Size >= 3)
    {
        int Rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

        if (Rank == 1)
        {
            double dVal = 7.7;
            MPI_Send(&dVal, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
        }
        if (Rank == 2)
        {
            int iVal = 3;
            MPI_Send(&iVal, 1, MPI_INTEGER, 0, 6, MPI_COMM_WORLD);
        }
        if (Rank == 0)
        {
            double dVal = 0;
            int iVal = 0;
            MPI_Status Status;

            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

            if (Status.MPI_TAG == 5)
            {
                MPI_Recv(&dVal, 1, MPI_DOUBLE, Status.MPI_SOURCE, 5, MPI_COMM_WORLD, &Status);
                std::cout << "Process 0 has received from process " << Status.MPI_SOURCE << ", dVal: " << dVal << std::endl;
            }
            if (Status.MPI_TAG == 6)
            {
                MPI_Recv(&iVal, 1, MPI_INTEGER, Status.MPI_SOURCE, 6, MPI_COMM_WORLD, &Status);
                std::cout << "Process 0 has received from process " << Status.MPI_SOURCE << ", iVal: " << iVal << std::endl;
            }

            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

            if (Status.MPI_TAG == 5)
            {
                MPI_Recv(&dVal, 1, MPI_DOUBLE, Status.MPI_SOURCE, 5, MPI_COMM_WORLD, &Status);
                std::cout << "Process 0 has received from process " << Status.MPI_SOURCE << ", dVal: " << dVal << std::endl;
            }
            if (Status.MPI_TAG == 6)
            {
                MPI_Recv(&iVal, 1, MPI_INTEGER, Status.MPI_SOURCE, 6, MPI_COMM_WORLD, &Status);
                std::cout << "Process 0 has received from process " << Status.MPI_SOURCE << ", iVal: " << iVal << std::endl;
            }
        }
    }
    else std::cout << "It's 3 or more processes need to be ran in parallel." << std::endl;

    MPI_Finalize();
}

void labwork_6_task_6()
{
    MPI_Init(NULL, NULL);

    int Size;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);

    if (Size >= 2)
    {
        int Rank;

        MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

        if (Rank == 0)
        {
            int N = 3 + (rand() % 10);
            double* A = new double[N];

            for (int k = 0; k < N; ++k)
                A[k] = k + 0.7;

            std::cout << "Process 0 sends the data: " << std::endl;

            for (int k = 0; k < N; ++k)
            {
                std::cout << A[k] << " ";
            }

            MPI_Send(A, N, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD);
            delete[] A;
        }
        if (Rank == 1)
        {
            MPI_Status Status;
            MPI_Probe(0, 5, MPI_COMM_WORLD, &Status);

            int N;
            MPI_Get_count(&Status, MPI_DOUBLE, &N);
            double* A = new double[N];
            MPI_Recv(A, N, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &Status);
            
            std::cout << "Process 1 received the data: " << std::endl;

            for (int k = 0; k < N; ++k)
                std::cout << A[k] << " ";

            delete[] A;
        }
    }
    else
    {
        std::cout << "2 processes are required to examine the performance." << std::endl;
    }

    MPI_Finalize();
}

void labwork_7_task_3()
{
    double a = 0, b = 0, c = 0;
    MPI_Request Request[4];
    MPI_Status Status[4];

    MPI_Init(NULL, NULL);

    int Size, Rank;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    int Posl = (Rank + 1) % Size;
    int Pred = Rank ? Rank - 1 : Size - 1;
    
    a = Rank + 0.7;
    
    MPI_Isend(&a, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &Request[0]);
    MPI_Irecv(&b, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &Request[1]);
    MPI_Isend(&a, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &Request[2]);
    MPI_Isend(&c, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &Request[3]);

    MPI_Waitall(4, Request, Status);

    std::cout << "Process " << Rank << " a=" << a << " b=" << b << " c=" << c << std::endl;

    MPI_Finalize();
}

void labwork_7_task_5()
{
    double a = 0, b = 0, c = 0;

    MPI_Request Request[4];
    MPI_Status Status[4];

    MPI_Init(NULL, NULL);

    int Size, Rank;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    int Posl = (Rank + 1) % Size;
    int Pred = Rank ? Rank - 1 : Size - 1;

    a = Rank + 0.7;

    MPI_Send_init(&a, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &Request[0]);
    MPI_Recv_init(&b, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &Request[1]);
    MPI_Send_init(&a, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &Request[2]);
    MPI_Recv_init(&c, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &Request[3]);

    for (int k = 1; k <= 3; ++k)
    {
        a += 0.1;
        MPI_Startall(4, Request);
        MPI_Waitall(4, Request, Status);
    }

    for (int k = 0; k < 4; ++k)
        MPI_Request_free(&Request[k]);
    std::cout << "Process " << Rank << " a=" << a << " b=" << b << " c=" << c << std::endl;

    MPI_Finalize();
}

bool IsFinished(bool* Prizn, int N)
{
    bool Res = true;

    for (int k = 0; k < N; ++k)
        Res = (Prizn[k]) && Res;

    return Res;
}

double Master(double* A, int N)
{
    double S = 0;

    for (int k = 0; k < N; ++k)
        S += A[k];

    return S;
}

void Slave(double* A, int N, int Rank)
{
    static int Count = 0;
    
    double* Tmp1 = new double[N];
    double* Tmp2 = new double[N];

    Tmp1[0] = cos(1.0 + Count + Rank);
    Tmp2[0] = cos(1.0 + Count - Rank);

    for (int k = 1; k < N; ++k)
    {
        Tmp1[k] = sin(Tmp1[k - 1] + Tmp2[k - 1]);
        Tmp2[k] = sin(Tmp1[k - 1] - Tmp2[k - 1]);
    }

    for (int k = 0; k < N; ++k)
    {
        A[k] = 0;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                A[k] += Tmp1[i] * Tmp2[i] / (1 + k + (i - j) * (i - j));
    }

    ++Count;
    delete[] Tmp2;
    delete[] Tmp1;
}

void labwork_7_task_6()
{
    MPI_Init(NULL, NULL);
    int Size, Rank;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    if (Size < 2)
        std::cout << "Minimum two processes required." << std::endl;
    else
    {
        if (Rank == 0)
        {
            bool* Prizn = new bool[Size - 1];
            double* A = new double[(Size - 1) * NN];
            double* R = new double[Size - 1];

            for (int k = 0; k < Size - 1; ++k)
                R[k] = 0;

            MPI_Request* Request = new MPI_Request[Size - 1];
            MPI_Status* Status = new MPI_Status[Size - 1];

            int* Indx = new int[Size - 1];
            double Tms = MPI_Wtime();

            for (int m = 0; m < NNN; ++m)
            {
                for (int k = 1; k < Size; ++k)
                {
                    Prizn[k - 1] = false;
                    MPI_Irecv(&A[NN * (k - 1)], NN, MPI_DOUBLE, k, 5, MPI_COMM_WORLD, &Request[k - 1]);
                }

                while (!IsFinished(Prizn, Size - 1))
                {
                    int Num;
                    MPI_Waitsome(Size - 1, Request, &Num, Indx, Status);

                    for (int k = 0; k < Num; ++k)
                    {
                        Prizn[Indx[k]] = true;
                        R[Indx[k]] += Master(&A[NN * Indx[k]], NN);
                    }
                }
            }

            Tms = MPI_Wtime() - Tms;

            for (int k = 0; k < Size - 1; ++k)
                std::cout << "R[" << k << "]=" << R[k] << std::endl;

            std::cout << "Time: " << Tms << " s" << std::endl;

            delete[] Indx;
            delete[] Status;
            delete[] Request;
            delete[] R;
            delete[] A;
            delete[] Prizn;
        }
        else
        {
            double* A = new double[NN];

            for (int m = 0; m < NNN; ++m)
            {
                Slave(A, NN, Rank);
                MPI_Send(A, NN, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
            }

            delete[] A;
        }
    }
    MPI_Finalize();
}

double MPI_Simps(double a, double b, int N, double Func(double), int& ProcID)
{
    double h = (b - a) / (2 * N);

    int k, NumProc;
    double S1, mpi_S1, S2, mpi_S2, Tmp;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &NumProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mpi_S1 = 0;
    mpi_S2 = 0;

    for (k = 1 + ProcID; k < N; k += NumProc)
    {
        Tmp = a + (2 * k - 1) * h;
        mpi_S1 += Func(Tmp);
        mpi_S2 += Func(Tmp + h);
    }

    MPI_Reduce(&mpi_S1, &S1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mpi_S2, &S2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    S1 += Func(b - h);

    return h * (Func(a) + Func(b) + 4.0 * S1 + 2.0 * S2) / 3.0;
}

#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

double MyFunction(double x)
{
    return 1.0 / (1.0 + x * x);
}

void labwork_8_task_3()
{
    int ProcID;
    double Tms = clock();
    double Intgr = MPI_Simps(0, 1000000, 1000000000, MyFunction, ProcID);

    Tms = (clock() - Tms) / CLOCKS_PER_SEC;

    if (!ProcID) {
        std::cout << "Time=" << Tms << " sec" << std::endl;
        std::cout.precision(8);
        std::cout << "Intgr=" << Intgr << ", M_PI=" << M_PI / 2.0 << std::endl;
    }
}

void labwork_8_task_5()
{
    MPI_Init(NULL, NULL);
    int Size;

    MPI_Comm_size(MPI_COMM_WORLD, &Size);

    int* Index = new int[Size];
    for (int k = 0; k < Size; ++k)
        Index[k] = Size - 1 + k;

    int* Edges = new int[2 * Size - 2];
    for (int k = 0; k < Size - 1; ++k)
    {
        Edges[k] = k + 1;
        Edges[Size - 1 + k] = 0;
    }

    MPI_Comm Gr_Comm;
    MPI_Graph_create(MPI_COMM_WORLD, Size, Index, Edges, true, &Gr_Comm);

    int Rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    int* Neighbors = new int[Size - 1];

    int Num;
    MPI_Graph_neighbors_count(Gr_Comm, Rank, &Num);
    MPI_Graph_neighbors(Gr_Comm, Rank, Num, Neighbors);

    for (int k = 0; k < Num; ++k)
    {
        int Rank1;
        MPI_Status Status;
        MPI_Sendrecv(&Rank, 1, MPI_INT, Neighbors[k], 5, &Rank1, 1, MPI_INT, Neighbors[k], 5, Gr_Comm, &Status);

        std::cout << "Process " << Rank << " co-operates with process " << Rank1 << std::endl;
    }

    MPI_Comm_free(&Gr_Comm);

    delete[] Neighbors;
    delete[] Edges;
    delete[] Index;

    MPI_Finalize();
}

struct TParticle
{
    double Coords[3];
    int Charge;
    double Mass;
    char Chr[8];
};

void labwork_8_task_7()
{
    TParticle Particle[3];

    MPI_Init(NULL, NULL);

    int Rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    MPI_Datatype TType[4] = { MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_CHAR };

    int BlockLen[4] = { 3, 1, 1, 8 };
    MPI_Aint Disp[4];
    MPI_Aint Base;

    MPI_Get_address(&Particle[0], &Base);
    MPI_Get_address(&(Particle->Coords), &Disp[0]);
    MPI_Get_address(&(Particle->Charge), &Disp[1]);
    MPI_Get_address(&(Particle->Mass), &Disp[2]);
    MPI_Get_address(&(Particle->Chr), &Disp[3]);

    for (int k = 0; k < 4; ++k)
        Disp[k] -= Base;

    MPI_Datatype MyParticle;
    MPI_Type_create_struct(4, BlockLen, Disp, TType, &MyParticle);
    MPI_Type_commit(&MyParticle);
    
    if (Rank == 0) {
        for (int k = 0; k < 3; k++)
            for (int j = 0; j < 8; j++)
                Particle[k].Chr[j] = 0;
        
        Particle[0].Coords[0] = 1.1;
        Particle[0].Coords[1] = 7.3;
        Particle[0].Coords[2] = 3.7;
        Particle[0].Charge = 1;
        Particle[0].Mass = 7.83;
        Particle[0].Chr[0] = 'P';
        Particle[0].Chr[1] = '1';
        Particle[1].Coords[0] = 1.1;
        Particle[1].Coords[1] = -7.3;
        Particle[1].Coords[2] = 3.7;
        Particle[1].Charge = -1;
        Particle[1].Mass = 1.83;
        Particle[1].Chr[0] = 'P';
        Particle[1].Chr[1] = '2';
        Particle[2].Coords[0] = -1.1;
        Particle[2].Coords[1] = 7.3;
        Particle[2].Coords[2] = -3.7;
        Particle[2].Charge = 0;
        Particle[2].Mass = 7.87;
        Particle[2].Chr[0] = 'P';
        Particle[2].Chr[1] = '3';
        
        std::cout << "Process 0 sends the data" << std::endl;
        for (int k = 0; k < 3; k++)
            std::cout << "Particle " << Particle[k].Chr << std::endl
            << "Coordinates " << Particle[k].Coords[0] << " "
            << Particle[k].Coords[1] << " " << Particle[k].Coords[2] << std::endl
            << "Mass " << Particle[k].Mass
            << " Charge " << Particle[k].Charge << std::endl;
        MPI_Send(&Particle[0], 3, MyParticle, 1, 5, MPI_COMM_WORLD);
    }
    if (Rank == 1) {
        MPI_Status St;
        MPI_Recv(&Particle[0], 3, MyParticle, 0, 5, MPI_COMM_WORLD, &St);

        std::cout << "Process 1 has accepted the data" << std::endl;

        for (int k = 0; k < 3; k++)
            std::cout << "Particle " << Particle[k].Chr << std::endl
            << "Coordinates " << Particle[k].Coords[0] << " "
            << Particle[k].Coords[1]
            << " " << Particle[k].Coords[2] << std::endl
            << "Mass " << Particle[k].Mass
            << " Charge " << Particle[k].Charge << std::endl;
    }

    MPI_Type_free(&MyParticle);
    MPI_Finalize();
}

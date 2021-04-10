# relief-mpi
Relief is a commonly used algorithm for feature selection in supervised machine learning. In this project, relief is implemented in MPI.
This is a course project of CMPE300( Analysis of Algorithms) at Boğaziçi University.

* Check [Relief Algorithm](https://medium.com/@yashdagli98/feature-selection-using-relief-algorithms-with-python-example-3c2006e18f83). (Medium)

## Execution
This project uses the Open MPI 4.0.3 implementation. For further information about Open MPI, please visit: https://www.open-mpi.org/ C++ is used as the programming language and the compiler is g++, version information: gcc version 9.3.0. 

_Compile:_ mpic++ -o mpi_relief ./mpi_relief.cpp

_Run:_ mpiexec --oversubscribe -np < P > mpi_relief < inputfile >

## Input Format
There are five parameters: P, N, A, M, T

P is the number of processors in the environment. N is number of instances in the input. A is number of features in the dataset. M is number of iterations that will be done in relief algorithm. T is number of features each processor will choose.


After these five parameters, there are N x (A+1) numbers in tabular format. The last number in each row corresponds to the label. 

An example input is present.

## Details
N instances are distributed to (P-1) slave processors. Each slave processor chooses T features and sends them to the master in a sorted format. The master merges features coming from all slave processors. Selection of instances that is sent to a slave processor is not random. Also, the slave processors don't choose random rows during M iterations, instead they work on rows ordered 1 to M.   

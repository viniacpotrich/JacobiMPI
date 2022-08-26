# JacobiMPI
Implementação do método de Jacobi com MPI
Eu optei por utilizar o MPI_Allgatherv, achei mais simples após testar com MPI_Allgather e MPI_Bcast.

#PARA COMPILAR
mpicc -o jacobi jacobiMPIV.c -lm

#PARA RODAR
//-np : é o número de threads
exemplo abaixo:
mpirun -np 1 ./jacobi 7 30 

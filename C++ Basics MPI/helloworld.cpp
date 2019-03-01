#include <iostream>
#include <mpi.h>

using namespace std;


int main(int argc, char ** argv){
int mynode, totalnodes;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
cout << "Hello world from process " << mynode;
cout << " of " << totalnodes << endl;
MPI_Finalize();
}

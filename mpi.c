#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include<sys/time.h>


const int MAX_STRING = 100;


int main(void) {
	struct timeval start, end; 
    gettimeofday(&start, NULL);
	char greeting[MAX_STRING];

	// no of processses
	int comm_sz;
	// my process rank
	int my_rank;

	MPI_Status status;
	MPI_Init(NULL,NULL);
	// & gives memory address of int variable
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	if (my_rank == 3) {
		sprintf(greeting,"Greetings from process %d of %d!\n", my_rank, comm_sz);
		MPI_Send(greeting, strlen(greeting)+1, MPI_CHAR,0,2,MPI_COMM_WORLD);
	}
	else if (my_rank == 2) {
		sprintf(greeting,"Greetings from process %d of %d!\n", my_rank, comm_sz);
		MPI_Send(greeting, strlen(greeting)+1, MPI_CHAR,0,2,MPI_COMM_WORLD);
	}
	else if (my_rank == 1) {
		int d = 9;
		int b = 4;
		int z = d*b;
		sprintf(greeting,"Greetings from process %d of %d!\n", my_rank, comm_sz);
		MPI_Send(greeting, strlen(greeting)+1, MPI_CHAR,0,2,MPI_COMM_WORLD);
	}
	else {
		printf(".........Greetings from process %d of %d!\n",my_rank,comm_sz);
		for (int q=1; q<comm_sz; q++) {
			MPI_Recv(greeting, MAX_STRING, MPI_CHAR, q,2, MPI_COMM_WORLD,&status);
		}
		printf("?????%s",greeting);

		// for (int q=1; q<comm_sz; q++) {
		// 	MPI_Recv(greeting, MAX_STRING, MPI_CHAR, MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
		//     printf("?????%s\n",greeting);
		// }
		
	}
	printf("FINALLY WE HAVE = %s\n",greeting);
	MPI_Finalize();

    gettimeofday(&end, NULL); 
    double time_taken;
    time_taken = (end.tv_sec - start.tv_sec) * 1e6; 
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6; 
    printf("Time Taken for process %d = %lf\n",my_rank, time_taken);


	return 0;
}
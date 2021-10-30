#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include<sys/time.h>

void Multiply_Serial(double *A, double *B, double *C, int m, int n, int p){
	int i, j, k;
	for (i = 0; i < m; i++){
		for (j = 0; j < p; j++){
			C[i*p + j] = 0;
			for (k = 0; k < n; k++)
				C[i*p + j] += A[i*n + k] * B[k*p + j];
		}
	}
}	

void Matrix_Multiply(double *A, double *B, double *C, int m, int n, int p, int comm_sz, int my_rank) {
	int i,j,k;
	int work = m/comm_sz;
	if (my_rank == comm_sz-1 && (my_rank+1)*work < m) {
		for (i=my_rank*work; i<m; i++) {
			for (j = 0; j < p; j++){
				C[i*p + j] = 0;
				for (k = 0; k < n; k++)
					C[i*p + j] += A[i*n + k] * B[k*p + j];
			}
		}
	}
	else {
		for (i=my_rank*work; i<(my_rank+1)*work; i++) {
			for (j = 0; j < p; j++){
				C[i*p + j] = 0;
				for (k = 0; k < n; k++)
					C[i*p + j] += A[i*n + k] * B[k*p + j];
			}
		}
	}
}

int IsEqual(double *A, double *B, int m, int n) {
		int i, j;
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			printf("(%0.2f,%0.2f),",A[i*m + n],B[i*m + n]);
			if (A[i*m + n] != B[i*m + n]) {
				printf("matrices not equal\n");
				return 0;
			}
		}
	}
	printf("matrices equal\n");
	return 1;
}

void print_matrix(double *A,int a) {
	// printf("size of this matrix = %d\n",sizeof(A[0]));
	printf("(");
	for (int i=0;i<a;i++) {
		printf("%0.6f, ",A[i]);
	}
	printf(")\n");

}

int main() {
	int N = 10000;
	int size;
	double* A = (double*) malloc(32*N*sizeof(double));
	double* B = (double*) malloc(32*N*sizeof(double));
    double* answer = (double*) malloc(N*N*sizeof(double));
	double* C = (double*) malloc(N*N*sizeof(double));
    double* answer_serial = (double*) malloc(N*N*sizeof(double));

    struct timeval start, end; 
	MPI_Init(NULL,NULL);
	int comm_sz;
	int my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// printf("bhag bsdk\n");
	if (my_rank != 0) {
		// printf("running process : %d\n", my_rank);
		MPI_Recv(&size,1,MPI_INT,0,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		// printf("size received = %d\n",size);
		MPI_Recv(A,32*N,MPI_DOUBLE,0,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(B,32*N,MPI_DOUBLE,0,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		int m=size, n = 32, p = size;
		for (int i=1;i<comm_sz;i++) {
			if (my_rank == i) {
				Matrix_Multiply(A,B,C,m,n,p,comm_sz,i);
				// printf("matrix C after process : %d\n",i),print_matrix(C,size*size);
				MPI_Send(C,size*size,MPI_DOUBLE,0,21,MPI_COMM_WORLD);	
			}
		}
	}

	else {
		printf("enter N = \n");
		scanf("%d",&size);
		double t0 = MPI_Wtime();
		for (int i=0; i<32*size; i++) {
			A[i] = (rand()/(double)RAND_MAX);
		}
		for (int i=0; i<32*size; i++) {
			B[i] = (rand()/(double)RAND_MAX);
		}
		// printf("matrices initiated\n");
		// printf("matrix A is :- \n"),print_matrix(A,32*size);
		// printf("matrix B is :- \n"),print_matrix(B,32*size);
		for (int i=1;i<comm_sz;i++) {
			MPI_Send(&size,1,MPI_INT,i,21,MPI_COMM_WORLD);
			// printf("sent size\n");
			MPI_Send(A,32*size,MPI_DOUBLE,i,21,MPI_COMM_WORLD);
			// printf("sent A\n");
			MPI_Send(B,32*size,MPI_DOUBLE,i,21,MPI_COMM_WORLD);
			// printf("sent B\n");	
		}
		int m=size, n = 32, p = size;
		Matrix_Multiply(A,B,answer,m,n,p,comm_sz,my_rank);
		// printf("matrix answer after process : %d\n",0),print_matrix(answer,size*size);
		for (int i=1;i<comm_sz;i++) {
			MPI_Recv(C,size*size,MPI_DOUBLE,i,21,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for (int j=0;j<size*size;j++) {
				if (C[j] != 0) {
					// printf("going to update answer matrix\n");
					answer[j] = C[j];
				}
			}
			// MPI_Recv(C,MAX_CAPACITY,MPI_double,i,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
		}

		// printf("answer using parallel is : "), print_matrix(answer,size*size);
		// gettimeofday(&start, NULL); 
		double t1 = MPI_Wtime();
		// Multiply_Serial(A,B,answer_serial,size,32,size);
		// double t2 = MPI_Wtime();
		// printf("time taken for serial = %f\n",t2-t1);
		printf("time taken for parallel = %f\n",t1-t0);
	} 

		
	MPI_Finalize(); 

	return 0;

}
#include "../include/block_matmul.h"

struct Config {
	/* MPI Files */
	MPI_File A_file, B_file, C_file;
	char *outfile;

	/* MPI Datatypes for matrix blocks */
	MPI_Datatype block;

	/* Matrix data */
	double *A, *A_tmp, *B, *C;

	/* Cart communicators */
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];
	int world_rank, world_size, grid_rank;
	int row_rank, row_size, col_rank, col_size;

	/* Full matrix dim */
	int A_dims[2];
	int B_dims[2];
	int C_dims[2];
	int matrix_size;

	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Config config;

void init_matmul(char *A_file, char *B_file, char *outfile)
{
	/* Copy output file name to configuration */
	MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);
	config.outfile = outfile;

	/* Get matrix size header */
	MPI_File_open(MPI_COMM_WORLD, A_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.A_file);
	MPI_File_open(MPI_COMM_WORLD, B_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.B_file);

	if(config.world_rank == 0){
		MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		config.matrix_size = config.A_dims[0]*config.A_dims[1];
	}

	MPI_Bcast(&config.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(config.world_rank != 0){
		config.A_dims[0] = sqrt(config.matrix_size);
		config.A_dims[1] = sqrt(config.matrix_size);
		config.B_dims[0] = sqrt(config.matrix_size);
		config.B_dims[1] = sqrt(config.matrix_size);
	}

	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
    MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Dims_create(config.world_size, 2, config.dim);

	/* Verify dim of A and B matches for matul and both are square*/
	if(config.A_dims[0] != config.A_dims[1] || config.B_dims[0] != config.B_dims[1] ||
		config.A_dims[1] != config.B_dims[0]) {
			printf("You have failed me for the last time young Timmy, wrong dimensions!");
			return;
	}

	config.local_dims[0] = config.A_dims[0] / config.dim[0];
	config.local_dims[1] = config.A_dims[1] / config.dim[1];
	config.local_size = config.local_dims[0] * config.local_dims[1];

	/* Create Cart communicator for NxN processes */
	int *period = (int*) malloc(2 * sizeof(int));
	period[0] = 1;
	period[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, period, 1, &config.grid_comm);

	/* Sub div cart communicator to N row communicator */
	int *rowRemains = (int*) malloc(2 * sizeof(int));

	rowRemains[0] = 0;
	rowRemains[1] = 1;//pure chansning att dom vill att row är true och col false
	MPI_Cart_sub(config.grid_comm, rowRemains, &config.row_comm);

	/* Sub div cart communicator to N col communicator */
	int *colRemains = (int*) malloc(2 * sizeof(int));
	colRemains[0] = 1;
	colRemains[1] = 0;//pure chansning att dom vill att row är false och col true
	MPI_Cart_sub(config.grid_comm, colRemains, &config.col_comm);

	MPI_Comm_rank(config.col_comm, &config.col_rank);
	MPI_Comm_rank(config.row_comm, &config.row_rank);
    MPI_Comm_size(config.col_comm, &config.col_size);
    MPI_Comm_size(config.row_comm, &config.row_size);

	/* Setup sizes of full matrices */

	/* Setup sizes of local matrix tiles */

	/* Create subarray datatype for local matrix tile */

	/* Create data array to load actual block matrix data */

	/* Set fileview of process to respective matrix block */

	/* Collective read blocks from files */

	MPI_Aint length = (sizeof(double))*config.local_size;
	MPI_Aint extent = config.world_size * length;
	MPI_Offset disp = (sizeof(int) * 2) + (config.world_rank * length);
	MPI_Datatype contig, filetype;
	MPI_Type_contiguous(config.local_size, MPI_DOUBLE, &contig);
	MPI_Type_create_resized(contig, 0, extent, &config.block);
	MPI_Type_commit(&config.block);
	MPI_File_set_view(config.A_file, disp, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	MPI_File_set_view(config.B_file, disp, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);

	double *testArray3 = (double*) malloc(config.local_size * sizeof(double));
	double *testArray4 = (double*) malloc(config.local_size * sizeof(double));
	double *testArray5 = (double*) malloc(config.local_size * sizeof(double));

	//double testArray3[config.local_size];
	//double testArray4[config.local_size];
	//double testArray5[config.local_size];

	


	MPI_File_read_all(config.A_file, testArray3, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read_all(config.B_file, testArray4, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

	config.A = testArray3;
	/*if(config.world_rank == (config.world_size - 1)){
	}*/
	//printf("C-ANSWERS %f\n", config.A[config.local_size-1]);
	//return;

	printf("%f\n\n", config.A[config.local_size-1]);
	config.B = testArray4;
	config.A_tmp = testArray5;

	// Initialize temp array
	for(int i = 0; i < config.local_size; i++){
		config.A_tmp[i] = config.A[i];
	}

	// Initialize C array to zero
	double *cArray = (double*) malloc(config.local_size * sizeof(double));
	for(int i = 0; i < config.local_size; i++){
		cArray[i] = 0.0;
	}
	config.C = cArray;

	/* Close data source files */
	MPI_File_close(&config.A_file);
	MPI_File_close(&config.B_file);
}

void cleanup_matmul()
{
	/* Rank zero writes header specifying dim of result matrix C */
	/*if(config.world_rank == (config.world_size - 1)){
		for(int i = 0; i < config.local_size; i++){
			printf("%f ", config.C[i]);
		}
	}*/
	printf("%f ", config.C[config.local_size-1]);

	/* Set fileview of process to respective matrix block with header offset */

	/* Collective write and close file */

	/* Cleanup */
}

void multiply_matrices()
{
  int i, j, k;
  // Loop nest optimized algorithm
  for (i = 0 ; i < config.local_dims[0] ; i++) {
    for (k = 0 ; k < config.local_dims[0] ; k++) {
      for (j = 0 ; j < config.local_dims[0] ; j++) {
        config.C[(i*config.local_dims[0]) + k] += config.A_tmp[(i*config.local_dims[0]) + j] * config.B[(j*config.local_dims[0]) + k];
		//printf("%f ",  config.C[(i*config.local_dims[0]) + k]);
	  }
    }
  }
}

void compute_fox()
{
	printf("%f\n\n", config.A[config.local_size-1]);
	
	/* Compute source and target for verticle shift of B blocks */
	int prev_rank, next_rank, recv_buffer[config.local_size];

	printf("Local size: %d\n", config.local_size);
	
    prev_rank = (config.col_rank + (config.col_size - 1)) % config.col_size;
    next_rank = (config.col_rank + 1) % config.col_size;
	for (int i = 0; i < config.dim[0]; i++) {
		printf("Loop %d\n", i);
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
		//MPI_Bcast(&config.A, 4, MPI_DOUBLE, 0, config.row_comm);
		MPI_Bcast(config.A_tmp, config.local_size, MPI_DOUBLE, (config.col_rank + i) % config.row_size, config.row_comm);

		//(diag_array[config.col_rank] + i) % config.row_size
		/* dgemm with blocks */
		multiply_matrices();

		//Reset A_tmp
		for(int j = 0; j < config.local_size; j++){
			config.A_tmp[i] = config.A[i];
		}
		
		/* Shfting block B upwards and receive from process below */
		if(config.col_rank == 0) {
        	MPI_Send(config.B, config.local_size, MPI_DOUBLE, prev_rank, 0, config.col_comm);
    	}
		double tempArray[config.local_size];
		MPI_Recv(tempArray, config.local_size, MPI_DOUBLE, next_rank, 0, config.col_comm, MPI_STATUS_IGNORE);
		//printf("Rank %d received %d from rank %d\n", rank, recv_rank, prev_rank);
		if(config.col_rank != 0) {
			//printf("Rank %d sending to%d\n", rank, next_rank);
			MPI_Send(config.B, config.local_size, MPI_DOUBLE, prev_rank, 0, config.col_comm);
		}
	}
	return;
}

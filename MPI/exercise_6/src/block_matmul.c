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
	int row_rank, row_size, col_rank, col_size, my_row, my_col;

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
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, period, 1, &(config.grid_comm));

	//int coordinates[2];
	int *coordinates = (int*) malloc(2 * sizeof(int));
	MPI_Comm_rank(config.grid_comm, &(config.grid_rank));
	MPI_Cart_coords(config.grid_comm, config.grid_rank, 2, coordinates);
  config.my_row = coordinates[0];
  config.my_col = coordinates[1];
	config.coords[0] = coordinates[0];
	config.coords[1] = coordinates[1];

	/* Sub div cart communicator to N row communicator */
	int *rowRemains = (int*) malloc(2 * sizeof(int));
	rowRemains[0] = 0;
	rowRemains[1] = 1;//pure chansning att dom vill att row är true och col false
	MPI_Cart_sub(config.grid_comm, rowRemains, &(config.row_comm));

	/* Sub div cart communicator to N col communicator */
	int *colRemains = (int*) malloc(2 * sizeof(int));
	colRemains[0] = 1;
	colRemains[1] = 0;//pure chansning att dom vill att row är false och col true
	MPI_Cart_sub(config.grid_comm, colRemains, &(config.col_comm));

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
	int starts[2] = {config.coords[0]*config.local_dims[0], config.coords[1]*config.local_dims[0]};
	int largeSize[2] = {config.A_dims[0], config.A_dims[0]};
	int smallSize[2] = {config.local_dims[0], config.local_dims[0]};
	MPI_Type_create_subarray(2, largeSize, smallSize, starts, MPI_ORDER_C, MPI_DOUBLE, &(config.block));
	MPI_Offset disp = (sizeof(int) * 2);

	/*Below commented code is for contiguous*/
	//MPI_Aint length = (sizeof(double))*config.local_dims[0];
	//MPI_Aint extent = config.dim[0] * length;// length;
	//MPI_Offset disp = (sizeof(int) * 2) + (sizeof(double))*((config.local_size*config.dim[0]*config.col_rank) + (config.local_dims[0] * config.row_rank));
	//MPI_Datatype contig;
	//MPI_Type_contiguous(config.local_dims[0], MPI_DOUBLE, &contig);
	//MPI_Type_create_resized(contig, 0, extent, &(config.block));

	MPI_Type_commit(&config.block);
	MPI_File_set_view(config.A_file, disp, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	MPI_File_set_view(config.B_file, disp, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);

	config.A = (double*) malloc(config.local_size * sizeof(double));
	config.B = (double*) malloc(config.local_size * sizeof(double));
	config.A_tmp = (double*) malloc(config.local_size * sizeof(double));

	MPI_File_read_all(config.A_file, config.A, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read_all(config.B_file, config.B, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

	// Initialize temp array
	for(int i = 0; i < config.local_size; i++){
		config.A_tmp[i] = config.A[i];
	}
	// Initialize C array to zero
	config.C = (double*) malloc(config.local_size * sizeof(double));
	for(int i = 0; i < config.local_size; i++){
		config.C[i] = 0.0;
	}

	/* Close data source files */
	MPI_File_close(&config.A_file);
	MPI_File_close(&config.B_file);
}

void cleanup_matmul()
{
	/* Rank zero writes header specifying dim of result matrix C */
	MPI_File_open(MPI_COMM_WORLD, "C_computed.answer", (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.C_file);
	if(config.world_rank == 0){
		MPI_File_write(config.C_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE); // Should be C_dims but the same
	}
	/* Set fileview of process to respective matrix block with header offset */
	MPI_Offset disp = (sizeof(int) * 2);
	//MPI_Offset disp = (sizeof(int) * 2) + (sizeof(double))*((config.local_size*config.dim[0]*config.col_rank) + (config.local_dims[0] * config.row_rank));
	MPI_File_set_view(config.C_file, disp, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	/* Collective write and close file */
	MPI_File_write(config.C_file, config.C, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&config.C_file);
	/* Cleanup */
}

void multiply_matrices(double* A, double* B, double* C) {
  int i, j, k;
  for (i = 0 ; i < config.local_dims[0] ; i++) {
    for (k = 0 ; k < config.local_dims[0] ; k++) {
      for (j = 0 ; j < config.local_dims[0] ; j++) {
        C[(i*config.local_dims[0]) + k] = C[(i*config.local_dims[0]) + k] + A[(i*config.local_dims[0]) + j] * B[(j*config.local_dims[0]) + k];
	  }
    }
  }
}

void compute_fox() {
	/* Compute source and target for verticle shift of B blocks */
	int source, dst, root;
	MPI_Cart_shift(config.col_comm, 0, -1, &source, &dst);
	for (int i = 0; i < config.dim[0]; i++) {
		root = (config.my_row + i) % config.dim[0];
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
		if(root == config.my_col) {
			MPI_Bcast(config.A, config.local_size, MPI_DOUBLE, root, config.row_comm);
			multiply_matrices(config.A, config.B, config.C);
		}
		else {
			MPI_Bcast(config.A_tmp, config.local_size, MPI_DOUBLE, root, config.row_comm);
			multiply_matrices(config.A_tmp, config.B, config.C);
		}
		/* dgemm with blocks */
		/* Shfting block B upwards and receive from process below */
		MPI_Sendrecv_replace(config.B, config.local_size, MPI_DOUBLE, dst, 0, source, 0, config.col_comm, MPI_STATUS_IGNORE);
	}
	return;
}

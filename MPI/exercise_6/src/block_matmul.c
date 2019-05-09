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
	
	printf("Rank %d; A_dims[0]: %d, A_dims[1]: %d\n", config.world_rank, config.A_dims[0], config.A_dims[1]);

	/* Get matrix size header */
	MPI_File_open(MPI_COMM_WORLD, A_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.A_file);
	MPI_File_open(MPI_COMM_WORLD, B_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.B_file);

	if(config.world_rank == 0){
	MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
	config.matrix_size = config.A_dims[0]*config.A_dims[1];
	}

	/* Broadcast global matrix sizes */
	MPI_Bcast(&config.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(config.world_rank != 0){
		config.A_dims[0] = sqrt(config.matrix_size);
		config.A_dims[1] = sqrt(config.matrix_size);
		config.B_dims[0] = sqrt(config.matrix_size);
		config.B_dims[1] = sqrt(config.matrix_size);
	}
	printf("Rank %d; Matrix size: %d\n", config.world_rank, config.matrix_size);
	printf("Rank %d; A_dims[0]: %d, A_dims[1]: %d\n", config.world_rank, config.A_dims[0], config.A_dims[1]);
	

	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
    MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Dims_create(config.world_size, 2, config.dim);

	printf("Rank %d; dim[0]: %d, dim[1]: %d\n", config.world_rank, config.dim[0], config.dim[1]);

	/* Verify dim of A and B matches for matul and both are square*/
	if(config.A_dims[0] != config.A_dims[1] || config.B_dims[0] != config.B_dims[1] ||
		config.A_dims[1] != config.B_dims[0]) {
			printf("You have failed me for the last time young Timmy, wrong dimensions!");
			return;
	}

	config.local_dims[0] = config.A_dims[0] / config.dim[0];
	config.local_dims[1] = config.A_dims[1] / config.dim[1];
	config.local_size = config.local_dims[0] * config.local_dims[1];
	printf("Rank %d; local_dims[0]: %d, local_dims[1]: %d, local_size: %d\n", config.world_rank, config.local_dims[0], config.local_dims[1], config.local_size);

	double testArrayTestie[100];
	if(config.world_rank == 1){
		printf("World rank: %d\n", config.world_rank);
		MPI_Offset temp = (2*sizeof(int)) + ((sizeof(double)*config.local_size) * config.world_rank);
		printf("Offset: %ld\n", (2*sizeof(int)) + ((sizeof(double)*config.local_size) * config.world_rank));
		//MPI_File_read_at(config.A_file, testArrayTestie, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
		MPI_File_read_at(config.A_file, temp, testArrayTestie, 100, MPI_DOUBLE, MPI_STATUS_IGNORE);
		for(int i = 0; i < 100; i++){
		//printf("config.A: %f\n", config.A[i]);
		printf("Element 1: %f\n", testArrayTestie[i]);
		}
	}
		

	/* Create Cart communicator for NxN processes */
	int period[2] = {1,1}; //chansning att det ska vara periodic topology
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, period, 1, &config.grid_comm);

	/* Sub div cart communicator to N row communicator */
	int rowRemains[2] = {0,1};//pure chansning att dom vill att row är true och col false
	MPI_Cart_sub(config.grid_comm, rowRemains, &config.row_comm);

	/* Sub div cart communicator to N col communicator */
	int colRemains[2] = {1,0};//pure chansning att dom vill att row är false och col true
	MPI_Cart_sub(config.grid_comm, colRemains, &config.col_comm);

	/* Setup sizes of full matrices */
	//chansar på att det är gjort

	/* Setup sizes of local matrix tiles */
	/* Outcommented for now
	config.local_dims[0] = config.A_dims[0] / config.dim[0];
	config.local_dims[1] = config.A_dims[1] / config.dim[1];
	config.local_size = config.local_dims[0] * config.local_dims[1];
	printf("Rank %d; local_dims[0]: %d, local_dims[1]: %d, local_size: %d\n", config.world_rank, config.local_dims[0], config.local_dims[1], config.local_size);
	*/
	/* Create subarray datatype for local matrix tile */
	MPI_Comm_rank(config.col_comm, &config.col_rank);
	MPI_Comm_rank(config.row_comm, &config.row_rank);
    MPI_Comm_size(config.col_comm, &config.col_size);
    MPI_Comm_size(config.row_comm, &config.row_size);



	int starts[2] = {config.col_rank, config.row_rank};
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, starts, MPI_ORDER_C, MPI_DOUBLE, &config.block); // MPI_INT? DOUBLE? WTF?
	MPI_Type_commit(&config.block);

	/* Create data array to load actual block matrix data */
	// Probably done

	

	/* Set fileview of process to respective matrix block */
	MPI_Offset viewOff = (sizeof(int) * 2) + (sizeof(double)*config.local_size * config.world_rank);
	MPI_File_set_view(config.A_file, viewOff, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	
	MPI_File_set_view(config.B_file, viewOff, MPI_DOUBLE, config.block, "native", MPI_INFO_NULL); // A file och B file?

	/* Collective read blocks from files */
	double testArray[config.local_size];
	double testArray2[config.local_size];
	MPI_File_read_all(config.A_file, testArray, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read_all(config.A_file, testArray2, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

	config.A = testArray;
	config.B = testArray2;

	for(int i = 0; i < config.local_size; i++){
		//printf("%f", testArray[i]);
	}
	//printf("\n");
	
	if(config.world_rank == 1){
		for(int i = 0; i < 100; i++){
		printf("config.A: %f\n", config.A[i]);
		}
	}
		
	

	/* Close data source files */
	MPI_File_close(&config.A_file);
	MPI_File_close(&config.B_file);
}

void cleanup_matmul()
{
	/* Rank zero writes header specifying dim of result matrix C */

	/* Set fileview of process to respective matrix block with header offset */

	/* Collective write and close file */

	/* Cleanup */
}

void compute_fox()
{

	/* Compute source and target for verticle shift of B blocks */

	for (int i = 0; i < config.dim[0]; i++) {
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */

		/* dgemm with blocks */

		/* Shfting block B upwards and receive from process below */

	}
}

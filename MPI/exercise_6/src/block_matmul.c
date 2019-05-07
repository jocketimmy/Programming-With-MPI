#include "block_matmul.h"

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
	config.outfile = outfile;

	/* Get matrix size header */
	MPI_File_open(MPI_COMM_WORLD, &A_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.A_file);
	MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_open(MPI_COMM_WORLD, &B_file, (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &config.B_file);
	MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
	config.matrix_size = config.A_dims[0]*config.A_dims[1];

	/* Broadcast global matrix sizes */
	MPI_Bcast(&config.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
    MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Dims_create(config.world_size, 2, &config.dim);


	/* Verify dim of A and B matches for matul and both are square*/
	if(config.A_dims[0] != config.A_dims[1] || config.B_dims[0] != config.B_dims[1] ||
		config.A_dims[1] != config.B_dims[0]) {
			printf("You have failed me for the last time young Timmy, wrong dimensions!");
			return;
	}

	/* Create Cart communicator for NxN processes */
	int period[2] = {1,1}; //chansning att det ska vara periodic topology
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, &period, 1, &grid_comm);

	/* Sub div cart communicator to N row communicator */
	int rowRemains[2] = {1,0};//pure chansning att dom vill att row är true och col false
	MPI_Cart_sub(config.grid_comm, rowRemains, &row_comm);

	/* Sub div cart communicator to N col communicator */
	int colRemains[2] = {0,1};//pure chansning att dom vill att row är false och col true
	MPI_Cart_sub(config.grid_comm, colRemains, &col_comm);

	/* Setup sizes of full matrices */
	//chansar på att det är gjort

	/* Setup sizes of local matrix tiles */
	config.local_dims[0] = config.A_dims[0] / config.dim[0];
	config.local_dims[1] = config.A_dims[1] / config.dim[1];

	/* Create subarray datatype for local matrix tile */
	int starts[2] = {0,0};
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, starts, MPI_ORDER_C, MPI_INT, config.block);

	/* Create data array to load actual block matrix data */

	/* Set fileview of process to respective matrix block */

	/* Collective read blocks from files */

	/* Close data source files */
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

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <mpi.h>

#include "pi.c"

#define NUM_ITER 100000000

void print_usage();

int main(int argc, char *argv[])
{
	int opt;
	int world_rank;
	int num_ranks;
	int count = 0, flip = 10000, seed = 1, result = 0;
	double pi = 0.0;
	char *filename = NULL;
	MPI_File fh;
	char resultString[4];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	int counts[num_ranks];
	MPI_Request request[num_ranks];
	MPI_Status status[num_ranks];
	flip = NUM_ITER/num_ranks;
	while ((opt = getopt(argc, argv, "s:f:o:")) != -1) {
		switch (opt) {
			case 's':
				seed = atoi(optarg);
				break;
			case 'f':
				flip = atoi(optarg);
				break;
			case 'o':
				filename = optarg;
				break;
			default:
				if (world_rank == 0) print_usage(argv[0]);
		}
	}

    double x, y, z;
    srand(seed*world_rank);
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < NUM_ITER/num_ranks; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));

        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count++;
        }
    }

	char out[100] = {world_rank + '0', ' '};
	//count/(NUM_ITER/num_ranks) + '0'
	//char worldrank[] = ;
	//char countre[] = {(count/(NUM_ITER/num_ranks)) + '0'};
	//strcat(out, worldrank);
	//strcat(out, countre);
	//char out[100] = {world_rank + '0'};
	//strcat(out, world_rank);
	//strcat(out, " ");
	//strcat(out, count/(NUM_ITER/num_ranks));
	//strcat(out, "\n");
	MPI_Offset MPIoff = 4 * world_rank;
	MPI_File_open(MPI_COMM_WORLD, "results.txt", (MPI_MODE_RDWR | MPI_MODE_CREATE), MPI_INFO_NULL, &fh);
	MPI_File_write_at(fh, MPIoff, &out, strlen(out), MPI_CHAR, MPI_STATUS_IGNORE);
	// Close the file
	MPI_File_close(&fh);

	
	MPI_Reduce(&count, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	printf("rank %d: %d / %d = %f\n", world_rank, count, flip, (double)count / (double)flip);
	if (world_rank == 0) {
		//count = 0;
		//for(int i = 0; i < num_ranks; i++) {
		//	count += counts[i];
		//}
		init_pi(seed*world_rank, filename);

		compute_pi(NUM_ITER, &result, &pi);

		//compute_pi(flip, &count, &pi);

		printf("pi: %f\n", pi);
	}

	cleanup_pi();
	MPI_Finalize();

	return 0;
}

void print_usage(char *program)
{
	fprintf(stderr, "Usage: %s [-s seed] [-f trials]\n", program);
	exit(1);
}
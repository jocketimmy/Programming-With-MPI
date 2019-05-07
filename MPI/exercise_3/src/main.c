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
	int count = 0, flip = 10000, seed = 1;
	double pi = 0.0;
	char *filename = NULL;


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
	MPI_Gather(&count, 1, MPI_INT, &counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	printf("rank %d: %d / %d = %f\n", world_rank, count, flip, (double)count / (double)flip);
	if (world_rank == 0) {
		count = 0;
		for(int i = 0; i < num_ranks; i++) {
			count += counts[i];
		}
		init_pi(seed*world_rank, filename);

		compute_pi(NUM_ITER, &count, &pi);

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

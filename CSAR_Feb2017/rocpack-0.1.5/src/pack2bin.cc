#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <err.h>
#include <errno.h>

#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <signal.h>
#include <sysexits.h>

#include "parse.h"
#include "particle.h"
#include "boundary.h"
#include "settings.h"
#include "gjk.h"

#include "sphere.h"

using namespace std;

/* global options */

int draw, draw_final, verbose;

/* global variables controlling time */

float t_start; /* time at start of simulation      */
float t_prev;  /* time of previous event           */
float t_curr;  /* time of current event            */
float t_next;  /* time of next event               */
float t_stop;  /* time when simulation should stop */

/* global variables controlling packing */

float target_fraction;

/* input/output variables */

char* input_fname;
char* output_fname;

/* command line options */

static struct option options[] = {
	{"help",      no_argument,           0, 'h'},
	{"version",   no_argument,           0, 'V'},
	{"verbose",   no_argument,    &verbose, 'v'},
	{"quiet",     no_argument,    &verbose, 'q'},
	{"layers",    required_argument,     0, 'l'},
	{0, 0, 0, 0}
};

void help()
{
	printf("Usage: pack [OPTIONS] INPUT_FILE OUTPUT_FILE\n\n");

	for (int i = 0; options[i].name != NULL; i++) {
		printf("\t-%c \t--%s\t%s\n\n",
			options[i].val, options[i].name,
			options[i].has_arg ? options[i].name : "");
	}
}

int main(int argc, char* argv[])
{
	if (argc == 1) { help(); exit(0); }

	/* global defaults */

	verbose = true;
	output_fname = NULL;
	target_fraction = 1.0;

	FILE* output;
	bool  periodic;  /* pack boundary is solid (false) or periodic (true) */

	int   nlayers = 1024;

	register_shape("sphere", new Sphere());

	parse_configuration_files();

	/* parse command line options */

	while (1) {
		int c, optidx = 0;

		c = getopt_long(argc, argv, "hl:qvV", options, &optidx);

		if (c == -1) break;

		switch (c) {

			case 0: /* long options */
				c = options[optidx].val;

				/* fall through */

			case 'q':
				verbose = false;
				break;

			case 'v':
				verbose = true;
				break;

			case 'V':
				printf("%s\n", PACKAGE_STRING);
				exit(0);

			case 'l':
				nlayers = (int) strtol(optarg, (char **)NULL, 10);

				if (nlayers < 0 || nlayers > 8192)
					errx(EX_DATAERR, "number of layers out of range\n");

				break;

			case 'h':
			case '?':
				help();
				exit(0);

			default:
				help();
				exit(1);
		}
	}

	switch (argc - optind) {
		case 0:
			errx(EX_NOINPUT, "no input file given.\n");

		case 1:
			errx(EX_NOINPUT, "no output file given.\n");

		case 2:
			input_fname  = argv[optind];
			output_fname = argv[optind+1];
			break;

		default:
			errx(EX_DATAERR, "too many arguments.\n");
	}

	parse_input_file(input_fname);

	Vector size = boundary->dimensions();

	/* create amira file */
	float pixel_size = size.max() / nlayers;

	int nx = size[0] / pixel_size;
	int ny = size[1] / pixel_size;
	int nz = size[2] / pixel_size;

	periodic = boundary->is_periodic();

	unsigned char* data = new unsigned char[nx * ny * nz];

	memset((void *)data, 0, sizeof(data)); /* make it all zeros */

	for (unsigned int n = 0; n < particle.size(); n++) {
		fprintf(stderr, "\r%5.1f%%", (100.0 * n) / particle.size());
		particle[n]->set_time(1.0);

		Point x = particle[n]->position();
		float r = particle[n]->bounding_radius();

		int imin = floor(((x[0] - r) + size[0]/2.0) / pixel_size);
		int jmin = floor(((x[1] - r) + size[1]/2.0) / pixel_size);
		int kmin = floor(((x[2] - r) + size[2]/2.0) / pixel_size);

		int imax =  ceil(((x[0] + r) + size[0]/2.0) / pixel_size);
		int jmax =  ceil(((x[1] + r) + size[1]/2.0) / pixel_size);
		int kmax =  ceil(((x[2] + r) + size[2]/2.0) / pixel_size);

		for (int k = kmin; k <= kmax; k++) {
			if (!periodic && (k < 0 || k >= nz)) continue;
			for (int j = jmin; j <= jmax; j++) {
				if (!periodic && (j < 0 || j >= ny)) continue;
				for (int i = imin; i <= imax; i++) {
					if (!periodic && (i < 0 || i >= nx)) continue;

					Point voxel((i+0.5)*pixel_size-size[0]/2.0,
								(j+0.5)*pixel_size-size[1]/2.0,
								(k+0.5)*pixel_size-size[2]/2.0);

					if (contains(*particle[n], voxel)) {
						if (periodic) {
							int ip = i, jp = j, kp = k;
							if (ip < 0) ip += nx; else if (ip >= nx) ip -= nx;
							if (jp < 0) jp += ny; else if (jp >= ny) jp -= ny;
							if (kp < 0) kp += nz; else if (kp >= nz) kp -= nz;
							data[(kp*ny + jp)*nx + ip] = particle[n]->tag();
						} else {
							data[(k*ny + j)*nx + i] = particle[n]->tag();
						}
					}
				}
			}
		}
	}

	fprintf(stderr, "\n");

	if (!(output = fopen(output_fname, "wb"))) {
		fprintf(stderr, "%s: %s\n", output_fname, strerror(errno));
		exit(1);
	}

	fwrite(&nx, sizeof(int), 1, output);
	fwrite(&ny, sizeof(int), 1, output);
	fwrite(&nz, sizeof(int), 1, output);
	fwrite(data, sizeof(unsigned char), nx*ny*nz, output);
	fclose(output);

	return 0;
}

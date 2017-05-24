/* RocStat - Statistics module for RocPack */

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "amira.h"
#include "random.h"

double        angle;
double        max_distance;
unsigned int  bins;
unsigned int  points_per_bin;
unsigned int  n_hotspots;
unsigned int  hotspot_radius;
unsigned char particle_bins;
char*         slice;
FILE *input,  *output;

void view();
void print_help();
void calculate_Sij();
void calculate_Sij_bins();
void calculate_Sijk();
void calculate_Sijk_angle();
void make_slice(char*);
void generate_hotspots(int, int);

int main(int argc, char* argv[]){

	int ch, command; char *slice = NULL;

	/* default options */
	command = 'i'; input  = stdin; output = stdout;

	hotspot_radius = 1;

	/* parse command line */
	while((ch = getopt(argc, argv, "23B:H:a:b:d:n:hio:pr:s:v")) != -1){
		switch(ch){
			case 's':
				slice = optarg;
			case 'v':
			case 'i':
			case 'p':
				command = ch;
				break;
			case '2':
				/* default values */
				command = ch; bins = 100; max_distance = 0.0; points_per_bin = 1000000; particle_bins = 0;
				break;
			case '3':
				/* default values */
				command = ch; bins = 100; max_distance = 0.0; points_per_bin =  200000; angle = 60 * M_PI / 180;
				break;
			case 'B': {
				int arg = (int) strtol(optarg, (char **)NULL, 10);

				if(arg < 0 || arg > 256){
					fprintf(stderr, "Error: %d: invalid number of particle bins, must be within 1..256.\n", arg);
				} else if(arg > 1)
					particle_bins = (unsigned char) arg;

				}
				break;
			case 'H':
				command = ch;
				n_hotspots = (int) strtol(optarg, (char **)NULL, 10);
				break;
			case 'a':
				angle = strtod(optarg, (char **)NULL);

				if(angle < 0.0 || angle > 180.0){
					command = 'a'; bins = 100; max_distance = 0.0;
					points_per_bin = 200000; angle = 0.0;
				}

				angle = angle * M_PI / 180.0; /* convert to radians */
				break;
			case 'b':
				bins = (int) strtol(optarg, (char **)NULL, 10);
				break;
			case 'd':
				max_distance = strtod(optarg, (char **)NULL);
				break;
			case 'n':
				points_per_bin = (int) strtol(optarg, (char **)NULL, 10);
				break;
			case 'o':
				if(!(output = fopen(optarg,"w"))){
					fprintf(stderr,"%s: %s: %s\n", argv[0], optarg, strerror(errno));
					exit(1);
				}
				break;
			case 'r':
				hotspot_radius = (unsigned int) strtol(optarg, (char **)NULL, 10);
				break;
			case 'h':
				print_help();
				exit(0);
				break;
			case '?':
			default:
				exit(1);
				break;
		}
	}

	if(argc - optind != 1){
		fprintf(stderr, "Error: No input files given.\n");
		exit(1);
	}

	load_amira_file(argv[optind], &mesh);

	/* set a reasonable default value if not set in command line */
	if(max_distance == 0.0)
		max_distance = 0.25*(mesh.size[0]+mesh.size[1]+mesh.size[2])/3.0;

	/* now do the useful stuff */
	switch(command){
		case 'i': {
			register int i; int nhits = 0;
			switch(mesh.data_size){
				case 1:
					for(i = 0; i < mesh.n_elements; i++)
						if(((unsigned char*)mesh.data)[i] > 0)
							nhits++;
					break;

				case 2:
					for(i = 0; i < mesh.n_elements; i++)
						if(((unsigned short*)mesh.data)[i] > 0)
							nhits++;
					break;
			}

			fprintf(stdout, "%s: %lux%lux%lu %s voxels (%.3fum resolution, %.6f volume fraction)\n",
				argv[optind], mesh.size[0], mesh.size[1], mesh.size[2],
				mesh.data_size == 1 ? "byte" : "ushort",
				mesh.psize[0], (double)(nhits)/mesh.n_elements);
		}
		break;

		case 'p': { /* volume fraction */
			register int i; int nhits = 0;
			switch(mesh.data_size){
				case 1:
					for(i = 0; i < mesh.n_elements; i++)
						if(((unsigned char*)mesh.data)[i] > 0)
							nhits++;
					break;

				case 2:
					for(i = 0; i < mesh.n_elements; i++)
						if(((unsigned short*)mesh.data)[i] > 0)
							nhits++;
					break;
				default:
					fprintf(stderr,"Error: unknown data type in Amira mesh file.\n");
					exit(1);
					break;
			}

			fprintf(stdout,"%.6f\n", (double)(nhits)/mesh.n_elements);
		}
			break;

		case '2': /* second order statistics */

			if(particle_bins == 0)
				calculate_Sij();
			else if(mesh.data_size > 1)
				calculate_Sij_bins();
			else {
				fprintf(stderr, "Warning: Input file contains single byte data, "
					" calculating Sij with no binning.\n");
				calculate_Sij();
			}

			break;

		case '3': /* third order statistics */

			calculate_Sijk();
			break;

		case 'a': /* third order statistics with varying angle */

			calculate_Sijk_angle();
			break;

		case 's': {/* write a slice file to output */
			make_slice(slice);
		} break;

		case 'v': {
			view(0);
		} break;

		case 'H':
			if (hotspot_radius == 0
					|| hotspot_radius > mesh.size[0]
					|| hotspot_radius > mesh.size[1]
					|| hotspot_radius > mesh.size[2]) {
				fprintf(stderr, "Error: bad hotspot radius\n");
				exit(1);
			}
			generate_hotspots(n_hotspots, hotspot_radius);
			break;

		default: /* unknown, this should never be reached... */
			fprintf(stderr,"%s: Unknown error.\n", argv[0]);
			exit(1);
			break;
	}

	free(mesh.data); fclose(input); fclose(output);

	return 0;
}


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <err.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <signal.h>

#include "config.h"
#include "parse.h"
#include "settings.h"
#include "packing.h"

#include "sphere.h"

using namespace std;

/* global options */

int verbose;
#ifdef HAVE_OPENGL
int draw, draw_final;
#endif

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
#ifdef HAVE_OPENGL
	{"draw",      no_argument,       &draw, 'd'},
	{"draw-final",no_argument, &draw_final, 'D'},
#endif
	{"fraction",  required_argument,     0, 'f'},
	{"levels",    required_argument,     0, 'l'},
	{"output",    required_argument,     0, 'o'},
	{"seed",      required_argument,     0, 's'},
	{0, 0, 0, 0}
};

void help()
{
	printf("Usage: pack [OPTIONS] INPUT_FILE\n\n");

	for (int i = 0; options[i].name != NULL; i++) {
		printf("\t-%c \t--%s\t%s\n\n",
			options[i].val, options[i].name,
			options[i].has_arg ? options[i].name : "");
	}
}

void signal_handler(int signo)
{
	if (signo == SIGINT || signo == SIGTERM) {
		target_fraction *= pow(t_curr/t_stop, 3.0);
		t_stop = t_curr;
	}
}

int main(int argc, char* argv[])
{
	if (argc == 1) { help(); exit(0); }

	signal(SIGINT, signal_handler);
	signal(SIGTERM, signal_handler);

	/* global defaults */

#ifdef HAVE_OPENGL
	draw = false;
	draw_final = false;
#endif
	verbose = true;
	output_fname = strdup("pack.out");
	target_fraction = 1.0;
	settings::set("growth_rate", 1.0);

	register_shape("sphere", new Sphere());

	parse_configuration_files();

	/* parse command line options */

	while (1) {
		int c, optidx = 0;

		c = getopt_long(argc, argv, "dDf:l:ho:s:qvV", options, &optidx);

		if (c == -1) break;

		switch (c) {

			case 0: /* long options */
				c = options[optidx].val;

				/* fall through */

#ifdef HAVE_OPENGL
			case 'd':
				draw = true;
				draw_final = true;
				break;
#endif

#ifdef HAVE_OPENGL
			case 'D':
				draw_final = true;
				break;
#endif

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
			{
				int levels = (int) strtol(optarg, NULL, 10);

				if (levels <= 0 || levels > 16)
					errx(1, "number of levels should be between 1 and 16: %d", levels);

				settings::set("hgrid_levels", levels);

				break;
			}

			case 'f':
				target_fraction = strtof(optarg, NULL);

				if (target_fraction <= 0.0 || target_fraction > 1.0)
					errx(1, "packing fraction out of range: %f", target_fraction);

				break;

			case 'o':
				output_fname = optarg;
				break;

			case 's':
				settings::set("cmdline_seed", (int) strtol(optarg, NULL, 10));
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

	input_fname = argv[optind];

	parse_input_file(input_fname);
	atexit(write_output_file);

	main_loop();

	return 0;
}

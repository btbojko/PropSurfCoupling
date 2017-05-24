#include "config.h"

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

#ifdef HAVE_GD
#include <gd.h>
#endif

#ifdef HAVE_TIFF
#include <tiffio.h>
#endif

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
	{"scale",     required_argument,     0, 'S'},
	{"slice",     required_argument,     0, 's'},
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

void write_txt(const char *fname, const char *data, int nx, int ny)
{
	FILE *output;

	if (!(output = fopen(fname, "w")))
		err(EX_OSERR, "%s", fname);

	for(int j = 0; j < ny; j++) {
		for(int i = 0; i < nx -1; i++)
			fprintf(output, "%c ", data[j*nx+i] == 1 ? '1' : '0');
		fprintf(output, "%c\n", data[j*nx + nx-1] == 1 ? '1' : '0');
	}

	fclose(output);
}

void write_amira(const char* fname, unsigned char *data, int nx, int ny, int nz, float scale)
{
    FILE* output;

	if (!(output = fopen(fname, "w")))
		err(EX_OSERR, "%s", fname);

    fprintf(output,"# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n\n");
    fprintf(output,"define Lattice %d %d %d\n\n", nx, ny, nz);
    fprintf(output,"Parameters {\n   Content \"%dx%dx%d byte, uniform coordinates\",\n", nx, ny, nz);
    fprintf(output,"   BoundingBox 0.000 %.3f 0.000 %.3f 0.000 %.3f\n   CoordType \"uniform\",\n}\n\n",
        scale * nx, scale * ny, scale * nz);
    fprintf(output,"Lattice { byte Data } @1\n\n@1\n");
    fwrite(data, sizeof(unsigned char), nx*ny*nz, output);
    fclose(output);
}

void write_png(const char *fname, const char *data, int nx, int ny)
{
#ifdef HAVE_GD
	FILE *output; gdImagePtr image = gdImageCreate(ny, nx);
	int black = gdImageColorAllocate(image,   0,   0,   0);
	int white = gdImageColorAllocate(image, 255, 255, 255);

	if (!(output = fopen(fname, "wb")))
		err(EX_OSERR, "%s", fname);

	for(int j = 0; j < ny; j++)
		for(int i = 0; i < nx -1; i++)
			gdImageSetPixel(image,  j,  i, data[j*nx+i] == 1 ? white : black);

	gdImagePng(image, output); gdImageDestroy(image); fclose(output);
#else
  errx(EINVAL, "Cannot write %s: no PNG support compiled in", fname);
#endif
}

void write_tiff(const char* fname, unsigned char *data, int nx, int ny, int nz, float resolution)
{
#ifdef HAVE_TIFF
    TIFF *output = TIFFOpen(fname, "w");

    if (!output) {
        fprintf (stderr, "Can't open tiff file for writing\n");
        exit(1);
    }

    for (int k = 0; k < nz; k++) {
        TIFFSetField(output, TIFFTAG_IMAGEWIDTH, nx);
        TIFFSetField(output, TIFFTAG_IMAGELENGTH, ny);
        TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(output, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        TIFFSetField(output, TIFFTAG_XRESOLUTION, resolution);
        TIFFSetField(output, TIFFTAG_YRESOLUTION, resolution);
        TIFFSetField(output, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER);
        TIFFSetField(output, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        TIFFSetField(output, TIFFTAG_PAGENUMBER, k, nz);

        for (int j = 0; j < ny; j++)
            TIFFWriteScanline(output, &data[(k*ny + j) * nx], j, 0);

        TIFFWriteDirectory(output);
    }

    TIFFClose(output);
#else
  errx(EINVAL, "Cannot write %s: no TIFF support compiled in", fname);
#endif
}

int main(int argc, char* argv[])
{
	if (argc == 1) { help(); exit(0); }

	/* global defaults */

	verbose = true;
	output_fname = NULL;
	target_fraction = 1.0;

	bool  periodic;  /* pack boundary is solid (false) or periodic (true) */

	int   nlayers = 1024;
	float scale   = 1.0; /* scaling factor S for the pack, i.e.  scale is 1:x */
	float cmd_scale = -1.0;
	char* slice_desc = NULL;

	register_shape("sphere", new Sphere());

	parse_configuration_files();

	/* parse command line options */

	while (1) {
		int c, optidx = 0;

		c = getopt_long(argc, argv, "hl:qs:S:vV", options, &optidx);

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

			case 's':
				slice_desc = optarg;
				break;

			case 'S':
				cmd_scale = strtof(optarg, (char **)NULL);

				if (cmd_scale < 0.0)
					errx(EX_DATAERR, "bad scaling factor.\n");

				break;

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

	float pixel_size = size.max() / nlayers;

	int nx = size[0] / pixel_size;
	int ny = size[1] / pixel_size;
	int nz = size[2] / pixel_size;

	periodic = boundary->is_periodic();

	if (slice_desc == NULL) {
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
								data[(kp*ny + jp)*nx + ip] = 255; // particle[i]->tag();
							} else {
								data[(k*ny + j)*nx + i] = 255; // particle[i]->tag();
							}
						}
					}
				}
			}
		}

		fprintf(stderr, "\n");

		/* adjust bounding box scaling */

		if (settings::isset("scaling"))
			scale = settings::get<float>("scaling");

		/* command line scaling has higher precedence */
		if (cmd_scale > 0.0)
			scale = cmd_scale;

		if (strcmp(strrchr(output_fname, '.'), ".am") == 0) {
			write_amira(output_fname, data, nx, ny, nz, scale * pixel_size);
		} else if (strcmp(strrchr(output_fname, '.'), ".tif")  == 0 ||
			     strcmp(strrchr(output_fname, '.'), ".tiff") == 0) {
			write_tiff(output_fname, data, nx, ny, nz, 10000.0 / (scale * pixel_size));
		} else {
			err(EX_DATAERR, "%s: unknown output file format (should be .am or .tiff)", output_fname);
		}

		delete[] data;
	} else {
		char* data;
		char* plane = strtok(slice_desc, ":=");
		char* value = strtok((char*)NULL,":=");

		switch (*plane) {
			case 'x':
			{
				int i = strtof(value, (char**)NULL) * nx;

				if (i < 0 || i >= nx)
					errx(EX_DATAERR, "slice location should be normalized to [0,1]\n");

				data = new char[nz*ny];

				for (unsigned int n = 0; n < particle.size(); n++) {
					particle[n]->set_time(1.0);
					Point x = particle[n]->position();
					float r = particle[n]->bounding_radius();

					int imin = floor(((x[0] - r) + size[0]/2.0) / pixel_size);
					int imax =  ceil(((x[0] + r) + size[0]/2.0) / pixel_size);

					int jmax =  ceil(((x[1] + r) + size[1]/2.0) / pixel_size);
					int jmin = floor(((x[1] - r) + size[1]/2.0) / pixel_size);

					int kmax =  ceil(((x[2] + r) + size[2]/2.0) / pixel_size);
					int kmin = floor(((x[2] - r) + size[2]/2.0) / pixel_size);

					/* particle does not intersect this slice */
					if (i < imin || i > imax) continue;

					for (int k = kmin; k <= kmax; k++) {
						if (!periodic && (k < 0 || k >= nz)) continue;
						for (int j = jmin; j <= jmax; j++) {
							if (!periodic && (j < 0 || j >= ny)) continue;

							Point voxel((i+0.5)*pixel_size-size[0]/2.0,
							            (j+0.5)*pixel_size-size[1]/2.0,
							            (k+0.5)*pixel_size-size[2]/2.0);

							if (contains(*particle[n], voxel)) {
								if (periodic) {
									int ip = i, jp = j, kp = k;
									if (ip < 0) ip += nx; else if (ip >= nx) ip -= nx;
									if (jp < 0) jp += ny; else if (jp >= ny) jp -= ny;
									if (kp < 0) kp += nz; else if (kp >= nz) kp -= nz;
									data[kp*nz+jp] = 1;
								} else {
									data[ k*nz+ j] = 1;
								}
							}
						}
					}
				}

				if (strcmp(strrchr(output_fname, '.'), ".png") == 0)
					write_png(output_fname, data, nz, ny);
				else
					write_txt(output_fname, data, nz, ny);

				break;
			}

			case 'y':
			{
				int j = strtof(value, (char**)NULL) * ny;

				if (j < 0 || j >= ny)
					errx(EX_DATAERR, "slice location should be normalized to [0,1]\n");

				data = new char[nz*nx];

				for (unsigned int n = 0; n < particle.size(); n++) {
					particle[n]->set_time(1.0);
					Point x = particle[n]->position();
					float r = particle[n]->bounding_radius();

					int imin = floor(((x[0] - r) + size[0]/2.0) / pixel_size);
					int imax =  ceil(((x[0] + r) + size[0]/2.0) / pixel_size);

					int jmax =  ceil(((x[1] + r) + size[1]/2.0) / pixel_size);
					int jmin = floor(((x[1] - r) + size[1]/2.0) / pixel_size);

					int kmax =  ceil(((x[2] + r) + size[2]/2.0) / pixel_size);
					int kmin = floor(((x[2] - r) + size[2]/2.0) / pixel_size);

					/* particle does not intersect this slice */
					if (j < jmin || j > jmax) continue;

					for (int k = kmin; k <= kmax; k++) {
						if (!periodic && (k < 0 || k >= nz)) continue;
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
									data[kp*nx+ip] = 1;
								} else {
									data[ k*nx+ i] = 1;
								}
							}
						}
					}
				}

				if (strcmp(strrchr(output_fname, '.'), ".png") == 0)
					write_png(output_fname, data, nx, nz);
				else
					write_txt(output_fname, data, nx, nz);

				break;
			}

			case 'z':
			{
				int k = strtof(value, (char**)NULL) * nz;

				if (k < 0 || k >= nz)
					errx(EX_DATAERR, "slice location should be normalized to [0,1]\n");

				data = new char[nx*ny];

				for (unsigned int n = 0; n < particle.size(); n++) {
					particle[n]->set_time(1.0);
					Point x = particle[n]->position();
					float r = particle[n]->bounding_radius();

					int imin = floor(((x[0] - r) + size[0]/2.0) / pixel_size);
					int imax =  ceil(((x[0] + r) + size[0]/2.0) / pixel_size);

					int jmax =  ceil(((x[1] + r) + size[1]/2.0) / pixel_size);
					int jmin = floor(((x[1] - r) + size[1]/2.0) / pixel_size);

					int kmax =  ceil(((x[2] + r) + size[2]/2.0) / pixel_size);
					int kmin = floor(((x[2] - r) + size[2]/2.0) / pixel_size);

					/* particle does not intersect this slice */
					if (k < kmin || k > kmax) continue;

					for (int j = jmin; j < jmax; j++) {
						if (!periodic && (j < 0 || j >= ny)) continue;
						for (int i = imin; i < imax; i++) {
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
									data[jp*nx+ip] = 1;
								} else {
									data[ j*nx+ i] = 1;
								}
							}
						}
					}
				}

				if (strcmp(strrchr(output_fname, '.'), ".png") == 0)
					write_png(output_fname, data, nx, ny);
				else
					write_txt(output_fname, data, nx, ny);

				break;
			}

			default:
				errx(EX_DATAERR, "plane must be either 'x', 'y', or 'z'\n");
		}

		delete data;
	}

	return 0;
}

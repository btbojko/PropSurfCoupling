/******************************************************************************/
/* rdf - calculate radial distribution function of point distributions        */
/*                                                                            */
/* Copyright 2014 Guilherme Amadio                                            */
/*                                                                            */
/* Permission is hereby granted, free of charge, to any person obtaining a    */
/* copy of this software to use, copy, modify, merge, publish, distribute,    */
/* sublicense, and/or sell copies of this software, and to permit persons     */
/* to whom the software is furnished to do so, subject to the following       */
/* conditions:                                                                */
/*                                                                            */
/* The above copyright notice and this permission notice shall be included in */
/* all copies or substantial portions of the software.                        */
/*                                                                            */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    */
/* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR       */
/* OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,      */
/* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR      */
/* OTHER DEALINGS IN THE SOFTWARE.                                            */
/*                                                                            */
/* To compile this software, use:                                             */
/*                                                                            */
/* g++ -O3 -msse -msse4 -fopenmp -o rdf rdf.cc                                */
/*                                                                            */
/* g++ 4.5 or more recent is recommended.                                     */
/*                                                                            */
/* The input file format is:                                                  */
/*                                                                            */
/* x1 y1 z1 r1 tag1 x2 y2 z2 r2 tag2 ... xn yn zn rn tagn                     */
/*                                                                            */
/* where tag is an integer. Particles with the same tag are grouped together  */
/* and gij can be calculated using the -t option. Run rdf --help or look at   */
/* the description below for more information.                                */
/******************************************************************************/

#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <err.h>
#include <errno.h>
#include <libgen.h>
#include <getopt.h>

#include <vector>

#ifdef HAVE_OPENMP
	#include <omp.h>
#endif

#if not defined(__SSE4_1__)
	#error "Cannot compile without support for SSE4.1 SIMD instructions"
#else
	#include <smmintrin.h> /* SSE4 intrinsics */
#endif

using namespace std;

typedef union {
	__m128 xyzr; float f[4];
} sse_vec;

inline __m128 _mm_abs_ps(__m128 v)
{
    return _mm_andnot_ps(_mm_set1_ps(-0.0f), v);
}

inline float _mm_norm3_ps(__m128 v)
{
	return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(v, v, 0x71)));
}

inline float get_component(__m128 m, const int i)
{
	sse_vec v; v.xyzr = m; return v.f[i];
}

void print_help(const char*);

int main(int argc, char* argv[]) {
	char *progname = basename(argv[0]);

	FILE *input  = NULL;
	FILE *output = stdout;

	unsigned long int *histogram;
	unsigned long int nbins = 512;

	bool scaled   = false; /* print r (false) or r/r_max (true) in output */
	bool verbose  = false; /* show progress on the terminal               */
	bool periodic = false; /* pack boundary is solid or periodic          */

	bool tagged = false;       /* for RDFs of different particles         */
	unsigned int tag1 = 0;     /* tag1 indicates the particle at center   */
	unsigned int tag2 = 0;     /* tag2 indicates neighbor particles       */
	vector<__m128*> particle;  /* vector to hold all particles (x,y,z,r)  */
	vector<__m128*> particle1; /* vector to hold particles tagged as tag1 */
	vector<__m128*> particle2; /* vector to hold particles tagged as tag2 */

	float R = 1.0, N;    /* maximum radius used in scaling and number density */
	float d_max = 0.0;   /* g(r) is calculated in the interval [0:d_max]      */
	float rho;           /* number density for neighbor particle type         */
	float W, H, D, L, V; /* width, height, depth, L = min(W, H, D), V = W*H*D */
	__m128 extent;       /* [W, H, D, 0.0] used for periodic RDFs             */

	if (argc < 2) {
		print_help(progname);
		exit(0);
	}

	{ /* parse command line arguments */

		int c, optind = 0;
		static struct option options[] = {
			{"output",		required_argument,	0, 'o'},
			{"bins",		required_argument,	0, 'b'},
			{"tag",			required_argument,	0, 't'},
			{"max-dist",	required_argument,	0, 'd'},
			{"periodic",	no_argument,		0, 'p'},
			{"scaled",		no_argument,		0, 's'},
			{"help",		no_argument,		0, 'h'},
			{"verbose",		no_argument,		0, 'v'},
			{0, 0, 0, 0},
		};

		while ((c = getopt_long(argc, argv,
			"b:d:ho:pst:v", options, &optind)) != -1) {

			switch (c) {

			case 'b':
				nbins = (unsigned long int) strtol(optarg, (char **)NULL, 10);
				break;

			case 'd':
				d_max = strtof(optarg, (char **)NULL);
				break;

			case 'p': periodic = true;
				break;

			case 's':
				scaled = true;
				break;

			case 't':
			{
				char* ctag1 = strsep(&optarg, ":");
				char* ctag2 = strsep(&optarg, ":");

				tagged = ctag1 && ctag2;

				if (ctag1)
					tag1 = (unsigned int) strtol(ctag1, (char **)NULL, 10);
				if (ctag2)
					tag2 = (unsigned int) strtol(ctag2, (char **)NULL, 10);

				break;
			}

			case 'v':
				verbose = true;
				break;

			case 'o':
				if (!(output = fopen(optarg, "w")))
					err(1, "%s", optarg);
				break;

			case 'h':
			case ':':
			case '?':
			default:
				print_help(progname);
				exit(0);
			}
		}
	}

	{ /* read in input file */

		if (argc - optind < 1) {
			input = stdin;
		} else {
			if ((input = fopen(argv[optind], "r")) == NULL)
					err(1, "%s", argv[optind]);
		}

		__m128 bbox_min = _mm_set1_ps(FLT_MAX);
		__m128 bbox_max = _mm_setzero_ps();

		while(!feof(input)) {
			float x, y, z, r;
			unsigned int c, tag;

			c = fscanf(input, "%f %f %f %f %u\n", &x, &y, &z, &r, &tag);

			if (c != 5)
				errx(1, "error reading input file");

			__m128 p = _mm_set_ps(r, z, y, x);

			bbox_min = _mm_min_ps(bbox_min, p);
			bbox_max = _mm_max_ps(bbox_max, p);

			if (!tagged) {
				if (tag == tag1)
					particle.push_back(new __m128(p));
			} else {
				if (tag == tag1)
					particle1.push_back(new __m128(p));

				if (tag == tag2)
					particle2.push_back(new __m128(p));
			}
		}

		fclose(input);

		if ((!tagged && particle.size() == 0) ||
		    ( tagged && (particle1.size() == 0 || particle2.size() == 0)))
			errx(1, "error: one or more tags do not exist in input file");

		/* calculate bounding box size */
		extent = _mm_sub_ps(bbox_max, bbox_min);

		/* round bounding box values from particles to compute volume */
		W = roundf(100.0*get_component(extent, 0)) / 100.0; /* width  */
		H = roundf(100.0*get_component(extent, 1)) / 100.0; /* height */
		D = roundf(100.0*get_component(extent, 2)) / 100.0; /* depth  */

		extent = _mm_set_ps(0.0, D, H, W);

		V = W * H * D;

		if (!tagged)
			rho = particle.size()  / V;
		else
			rho = particle2.size() / V;

		L = W < D ? (W < H ? W : H) : (H < D ? H : D); /* min(W, H, D) */

		if (periodic) { /* check if centered at the origin */
			__m128 center = _mm_add_ps(bbox_max, bbox_min);

			if (_mm_norm3_ps(center) > L/10.0)
				errx(1, "error: periodic packings "
					"must be centered at the origin");
		}

		if (scaled) {
			unsigned int i;

			if (!tagged) {
				R = get_component(*particle[0], 3);
				for (i = 1; i < particle.size(); i++) {
					float r = get_component(*particle[i], 3);
					if (r > R) R = r;
				}
			} else {
				float r1 = get_component(*particle1[0], 3);
				float r2 = get_component(*particle2[0], 3);

				for (i = 1; i < particle1.size(); i++) {
					float r = get_component(*particle1[i], 3);
					if (r > r1) r1 = r;
				}

				for (i = 1; i < particle2.size(); i++) {
					float r = get_component(*particle2[i], 3);
					if (r > r2) r2 = r;
				}

				R  = r1 > r2 ? r1 : r2;
			}

			d_max *= R;
		}

		if (d_max == 0.0) { /* unset */
			d_max = periodic ? L / 3.0 : L / 4.0;
		} else if (d_max > (periodic ? L / 3.0 : L / 4.0)) {
			warnx("warning: maximum distance too large, decreasing to %s.",
				periodic ? "L/3" : "L/4");
			d_max = periodic ? L / 3.0 : L / 4.0;
		}
	}

	{ /* calculate the radial distribution function */

		unsigned int i, j;
		unsigned int nskip = 0;
		unsigned int progress = 0;
		float dr = d_max / nbins;
		float inv_dr = 1.0 / dr;

		histogram = new unsigned long int[nbins];
		bzero((void *)histogram, nbins * sizeof(unsigned long int));

		if (periodic) {
			if (!tagged) {
				#pragma omp parallel for default(shared) private(i,j)
				for (i = 0; i < particle.size(); i++) {
					__m128 p = *particle[i];

					if (verbose && progress++ % 100)
						fprintf(stderr, "\r%5.1f%%",
							100.0 * progress / particle.size());

					for (j = i + 1; j < particle.size(); j++) {
						/* compute distance between minimum images */
						__m128 tmp1 = _mm_abs_ps(_mm_sub_ps(p, *particle[j]));
						__m128 tmp2 = _mm_abs_ps(_mm_sub_ps(tmp1, extent));
						__m128 diff = _mm_min_ps(tmp1, tmp2);
						float  dist = _mm_norm3_ps(diff);
						unsigned int bin = (unsigned int) floor(dist * inv_dr);
						if (bin < nbins)
							#pragma omp atomic
							histogram[bin] += 2; /* count ij and ji */
					}
				}
			} else { /* tagged, periodic */
				#pragma omp parallel for default(shared) private(i,j)
				for (i = 0; i < particle1.size(); i++) {
					__m128 p = *particle1[i];

					if (verbose && progress++ % 100)
						fprintf(stderr, "\r%5.1f%%",
							100.0 * progress / particle1.size());

					for (j = 0; j < particle2.size(); j++) {
						/* compute distance between minimum images */
						__m128 tmp1 = _mm_abs_ps(_mm_sub_ps(p, *particle2[j]));
						__m128 tmp2 = _mm_abs_ps(_mm_sub_ps(tmp1, extent));
						__m128 diff = _mm_min_ps(tmp1, tmp2);
						float  dist = _mm_norm3_ps(diff);
						unsigned int bin = (unsigned int) floor(dist * inv_dr);
						if (bin < nbins)
							#pragma omp atomic
							histogram[bin]++;
					}
				}

			}
		} else { /* not tagged, not periodic */
			if (!tagged) {
				#pragma omp parallel for default(shared) private(i,j)
				for (i = 0; i < particle.size(); i++) {
					__m128 p = *particle[i];

					if (verbose)
						fprintf(stderr, "\r%5.1f%%",
							100.0 * progress++ / particle.size());

					/* skip particles too close to the boundary */
					if (_mm_norm3_ps(p) + d_max > L / 2.0) {
						nskip++;
						continue;
					}

					for (j = i+1; j < particle.size(); j++) {
						__m128 diff = _mm_sub_ps(p, *particle[j]);
						float  dist = _mm_norm3_ps(diff);
						unsigned int bin = (unsigned int) floor(dist * inv_dr);
						if (bin < nbins)
							#pragma omp atomic
							histogram[bin] += 2; /* count ij and ji */
					}
				}
			} else { /* tagged, not periodic */
				#pragma omp parallel for default(shared) private(i,j)
				for (i = 0; i < particle1.size(); i++) {
					__m128 p = *particle1[i];

					if (verbose)
						fprintf(stderr, "\r%5.1f%%",
							100.0 * progress++ / particle1.size());

					/* skip particles too close to the boundary */
					if (_mm_norm3_ps(p) + d_max > L / 2.0) {
						nskip++; continue;
					}

					for (j = 0; j < particle2.size(); j++) {
						__m128 diff = _mm_sub_ps(p, *particle2[j]);
						float  dist = _mm_norm3_ps(diff);
						unsigned int bin = (unsigned int) floor(dist * inv_dr);
						if (bin < nbins)
							#pragma omp atomic
							histogram[bin]++;
					}
				}
			}
		}

		if (verbose)
			fprintf(stderr, "\n");

	/* write output file */

		if (!tagged) {
			N = rho * (particle.size()  - nskip);
		} else {
			N = rho * (particle1.size() - nskip);
		}

		for (unsigned int i = 1; i < nbins; i++) {
			float r = (i + 0.5) * dr / R;
			float v = 4.0/3.0 * M_PI *
				(pow((float)(i+1), 3) - pow((float)(i), 3)) * pow(dr, 3);
			fprintf(output, "%lf %lf\n", r, ((float)(histogram[i])) / (N * v));
		}

		delete histogram; fclose(output);
	}


	return 0;
}

void print_help(const char* progname) {

static const char* help =
	"  Options:\n\n"
	"\t-b --bins\tBINS\tset number of bins to BINS (default: 512)\n\n"
	"\t-d --max-dist\tDIST\tcalculate g(r) in the interval [0:DIST]\n"
	"\t             \t    \t(default: L/3 for periodic, L/4 otherwise)\n\n"
	"\t-h --help\t\tprint help and exit\n\n"
	"\t-o --output\tFILE\tsave output to FILE (default: stdout)\n\n"
	"\t-p --periodic\t\tuse algorithm for a periodic pack (default: false)\n\n"
	"\t-s --scaled\t\toutput r as r/R, where R = max(ri, rj) (default: false)\n\n"
	"\t-t --tag\ti[:j]\tcalculate g_ii(r) or g_ij(r) if j is set (default: 0)\n\n"
	"\t-v --verbose\t\tshow progress meter of calculation (default: false)\n\n";

	printf("Usage: %s [options] input_file\n\n%s", progname, help);
}

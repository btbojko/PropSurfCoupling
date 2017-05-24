#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <cstring>

#include <err.h>
#include <time.h>
#include <sysexits.h>
#include <sys/types.h>
#include <dirent.h>

#include "token.h"
#include "random.h"

#include "shape.h"
#include "cylinder.h"
#include "ellipsoid.h"
#include "polyhedron.h"

#include "particle.h"
#include "boundary.h"
#include "settings.h"
#include "collision.h"

static token_t ctoken;

static int line;
static int header;
static FILE* input;
extern char* input_fname;
extern char* output_fname;
extern float t_curr;

void syntax_error(const char* msg)
{
    errx(EX_DATAERR, "syntax error: line %d: %s\n", line+1, msg);
}

void shift()
{
	ctoken = next_token(input, line);
}

int accept(int token)
{
    if (ctoken.type == token) {
        shift();
		return true;
    }

    return false;
}

void expect(int token)
{
    if (accept(token) == false) {
		char expected[64], found[64];

		strncpy(expected, token_description(token), sizeof(expected));
		strncpy(found, token_description(ctoken.type), sizeof(found));

		errx(EX_DATAERR, "syntax error: line %d:"
			" unexpected token: \"%s\" (expected \"%s\")\n",
			line, found, expected);
	}
}

int get_integer()
{
	if (ctoken.type != NUMBER)
		syntax_error("expected an integer");

	int result = (int) ctoken.value.fval;

	shift();

	return result;
}

float get_float()
{
	if (ctoken.type != NUMBER)
		syntax_error("expected a floating point number");

	float result = ctoken.value.fval;

	shift();

	return result;
}

Vector get_vector()
{
	float x, y, z;

	expect('<');

	x = get_float();
	accept(',');
	y = get_float();
	accept(',');
	z = get_float();

	expect('>');

	return Vector(x, y, z);
}

char* get_identifier()
{
	if (ctoken.type != IDENTIFIER)
		syntax_error("expected an identifier");

	assert(ctoken.value.str != NULL);

	char *identifier = strdup(ctoken.value.str);

	shift();

	return identifier;
}

/*** grammar rules start here ***/

/*** settings ***/

void parse_settings(void)
{
	while (ctoken.type == SET || ctoken.type == UNSET) {
		while (accept(SET)) {
			char *id = get_identifier();

			accept('=');

			switch (ctoken.type) {

				case ';':
					settings::set(id, true);
					break;

				case NUMBER:
					settings::set(id, ctoken.value.fval);
					break;

				case STRING:
				case IDENTIFIER:
					settings::set(id, ctoken.value.str);
					break;

				default:
					syntax_error("expected a number or a string");
					break;
			}

			free(id);

			shift();

			expect(';');
		}

		while (accept(UNSET)) {
			char *id = get_identifier();

			settings::unset(id);

			free(id);

			expect(';');
		}
	}
}

/*** new shape definition ***/

void parse_new_shape(void)
{
	const char *name = get_identifier();

	expect('{');

		if (strcmp(ctoken.value.str, "ellipsoid") == 0) {

			shift();

			expect('{');
				float a, b, c;
				a = get_float();
				accept(',');
				b = get_float();
				accept(',');
				c = get_float();
			expect('}');

			register_shape(name, new Ellipsoid(name, a, b, c));

			expect('}'); /* close shape definition */

		} else if (strcmp(ctoken.value.str, "cylinder") == 0) {

			shift();

			expect('{');
				float r, h;
				r = get_float();
				accept(',');
				h = get_float();
			expect('}');

			register_shape(name, new Cylinder(name, r, h));

			expect('}'); /* close shape definition */

		} else {
			std::vector<Point> vertices;
			std::vector<std::vector<int> > faces;

			if (strcmp(ctoken.value.str, "vertices") != 0)
				syntax_error("error in shape definition");
			else
				shift();

			expect('{');

				while (ctoken.type == '<') {
					vertices.push_back(get_vector());
					accept(',');
				}

			expect('}');

			if (strcmp(ctoken.value.str, "faces") != 0)
				syntax_error("error in shape definition");
			else
				shift();

			expect('{');

				while (ctoken.type == '[') {
					std::vector<int> face;

					shift();

					while (ctoken.type == NUMBER) {
						assert(ctoken.value.fval < vertices.size());
						face.push_back(get_integer());
						accept(',');
					}

					expect(']');
					accept(',');

					faces.push_back(face);
				}

			expect('}');

		expect('}');

		register_shape(name, new Polyhedron(name, vertices, faces));
	}
}

void parse_boundary(void)
{
	if (!accept(BOUNDARY))
		syntax_error("no boundary definition");

	expect('{');

	if (ctoken.type != IDENTIFIER)
		syntax_error("unexpected token, expected boundary type");

	char* boundary_name = get_identifier();

	if (strcmp(boundary_name, "box") == 0) {
		bool periodic = false;
		float width  = get_float();
		float height = get_float();
		float depth  = get_float();

		if (ctoken.type == IDENTIFIER) {
			if (strcmp(ctoken.value.str, "periodic") == 0) {
				periodic = true;
			} else if (strcmp(ctoken.value.str, "solid") == 0) {
				periodic = false;
			} else {
				syntax_error("expected 'solid' or 'periodic'");
			}

			shift();
		}

		if (periodic) {
			period = Vector(width, height, depth);
			boundary = new PeriodicBoundary(width, height, depth);
		} else {
			boundary = new BoxBoundary(width, height, depth);
		}
#ifndef EDMC
	} else if (strcmp(boundary_name, "sphere") == 0) {
		float radius = get_float();
		boundary = new SphericBoundary(radius);
	} else if (strcmp(boundary_name, "torus") == 0) {
		float R = get_float();
		float r = get_float();
		boundary = new TorusBoundary(R, r);
	} else if (strcmp(boundary_name, "cylinder") == 0) {
		float r = get_float();
		float h = get_float();
		boundary = new CylindricBoundary(r, h);
	} else if (strcmp(boundary_name, "annulus") == 0) {
		float r = get_float();
		float R = get_float();
		float h = get_float();
		boundary = new AnnulusBoundary(r, R, h);
	} else if (strcmp(boundary_name, "annulus_wedge") == 0) {
		float r = get_float();
		float R = get_float();
		float h = get_float();
		float t = get_float();
		if (t < 0.0 || t > 180.0)
			syntax_error("wedge angle must be between 0.0 and 180.0 degrees");
		boundary = new AnnulusWedgeBoundary(r, R, h, t*M_PI/180.0);
#endif
	} else {
		char msg[64];
		sprintf(msg, "unknown boundary type: %s", boundary_name);
		syntax_error(msg);
	}

	expect('}');
}

void parse_init_bare(void)
{
	int n, rng_seed;
	float max_size = 0.0;
	float min_size = INFINITY;
	Shape *shape; char *shape_name;
	rng_distribution* sdist;

	/* default size distribution */
	sdist = new constant_dist(1.0);

	/* if seed is set in input file, use it */
	if (settings::isset("seed"))
		rng_seed = settings::get<int>("seed");
	else
		rng_seed = time(NULL);

	/* command-line seed overwrites input file seed */
	if (settings::isset("cmdline_seed"))
		rng_seed = settings::get<int>("cmdline_seed");

	srand48(rng_seed);
	sdist->seed(rng_seed);

	while (accept(CREATE)) {
		n = get_integer();
		shape_name = get_identifier();

		unsigned int tag = 0;
		float r = 0.1529, g = 0.4911;
		float b = 0.7843, a = 1.0000;

		while (true) {
			if (accept(SIZE)) {
				if (sdist) delete sdist;
				switch (ctoken.type) {
					case NUMBER: {
						float s = get_float();

						if (s < 0.0)
							syntax_error("size must be a positive number");

						sdist = new constant_dist(s);
						break;
					}
					case UNIFORM:
					{
						shift();
						float a = get_float();
						float b = get_float();

						if (a < 0.0 || b < 0.0 || a > b)
							errx(EX_DATAERR, "syntax error: line %d: "
								"bad size interval: [%f:%f]\n", line+1, a, b);

						sdist = new uniform_dist(a, b);
						break;
					}
					case LOGNORMAL:
					{
						shift();
						float mean  = get_float();
						float sigma = get_float();

						if (sigma < 0.0)
							syntax_error("standard deviation be positive");

						sdist = new lognormal_dist(mean, sigma);
						break;
					}
					case GAUSSIAN:
					{
						shift();
						float mean  = get_float();
						float sigma = get_float();

						if (mean  < 0.0 || sigma < 0.0)
							syntax_error("mean size and standard deviation"
								" must be both positive");

						sdist = new normal_dist(mean, sigma);
						break;
					}
#if 0
					case GAMMA:
					{
						shift();
						float a = get_float();
						float b = get_float();

						if (a < 0.0 || b < 0.0)
							syntax_error("alpha and beta must be both positive");

						sdist = new gamma_dist(a, b);
						break;
					}
#endif
					case WEIBULL:
					{
						shift();
						float a = get_float();
						float b = get_float();

						if (a < 0.0 || b < 0.0)
							syntax_error("a and b must be both positive");

						sdist = new weibull_dist(a, b);
						break;
					}
					case CAUCHY:
					{
						shift();
						float a = get_float();
						float b = get_float();

						if (a < 0.0 || b < 0.0)
							syntax_error("a and b must be both positive");

						sdist = new cauchy_dist(a, b);
						break;
					}
				}
			} else if (accept(TAG)) {
				tag = get_integer();
			} else if (accept(COLOR)) {
				switch (ctoken.type) {
					case NUMBER:
						r = get_float();
						g = get_float();
						b = get_float();
						break;

					case IDENTIFIER:
						shift();
						r = drand48();
						g = drand48();
						b = drand48();
						break;
				}

			} else {
				break;
			}
		}

		expect(';');

		shape = create_shape(shape_name);

		assert(shape != NULL);

		for (int i = 0; i < n; i++) {
			Particle *p = new Particle(shape);
			float size = (*sdist)();
			set_min(min_size, size);
			set_max(max_size, size);
			p->set_growth_rate(size);
			p->set_tag(tag);
			p->set_color(r, g, b, a);
			particle.push_back(p);
		}

		free(shape_name);
	}

	settings::set("min_size", min_size);
	settings::set("max_size", max_size);
}

void parse_init_pack(void)
{
	float max_size = 0.0;
	float min_size = INFINITY;

	/* empty for now */
	while(ctoken.type == IDENTIFIER) {
		char* shape_name = get_identifier();
		Shape *shape = create_shape(shape_name);

		Particle *p = new Particle(shape);

		expect('{');

		while (true) {
			if (accept(TRANSLATE)) {
				p->translate(get_vector());
			} else if (accept(ROTATE)) {
				Vector axis = get_vector();
				float angle = get_float();
				p->rotate(Quaternion(axis, rad(angle)));
			} else if (accept(SIZE) || accept(SCALE)) {
				float size = get_float();
				set_min(min_size, size);
				set_max(max_size, size);
				p->set_growth_rate(size);
			} else if (accept(TAG)) {
				p->set_tag(get_integer());
			} else if (accept(COLOR)) {
				float r = 0.0f, g = 0.0f, b = 0.0f, a = 1.0f;
				switch (ctoken.type) {
					case NUMBER:
						r = get_float();
						g = get_float();
						b = get_float();
						break;

					case IDENTIFIER:
						shift();
						r = drand48();
						g = drand48();
						b = drand48();
						break;
				}
				p->set_color(r, g, b, a);
			} else {
				break;
			}
		}

		expect('}');

		particle.push_back(p);
	}

	settings::set("min_size", min_size);
	settings::set("max_size", max_size);
	settings::set("noinit", true);
}

void parse(const char* filename);
void parse_config(const char* filename);

void import(const char* filename)
{
	static int depth = 0;
	/* save current file and position */

	if (depth == 0) {
		FILE* cur_input = input;
		int cur_line = line;

		depth++;

		parse_config(filename);

		depth--;

		/* restore current file and position */

		input = cur_input;
		line = cur_line;

		shift(); /* shift previous EOF token */
	} else {
		syntax_error("recursive imports not supported");
	}
}

/* configuration files cannot define boundary and particles */

void parse_config(const char* filename)
{
	/* open input file */
	if ((input = fopen(filename, "r")) == NULL)
		return;

	line = 0; shift(); /* load first token */

	parse_settings();

	while (accept(SHAPE))
		parse_new_shape();

	parse_settings();

	fclose(input);
}

void parse_configuration_files()
{
	char path[512], fullpath[512], *home = getenv("HOME");

	DIR *dir; struct dirent *ent;
	/* scan system directory for shape files */
	if ((dir = opendir("/usr/share/rocpack/shapes")) != NULL) {
		while ((ent = readdir (dir)) != NULL) {
			if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, ".."))
				continue;
			sprintf(fullpath, "%s/%s",
				"/usr/share/rocpack/shapes", ent->d_name);
			parse_config(fullpath);
		}
		closedir(dir);
	}

	/* scan home directory for shape files */
	sprintf(path, "%s/.rocpack/shapes", home);
	if ((dir = opendir(path)) != NULL) {
		while ((ent = readdir (dir)) != NULL) {
			if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, ".."))
				continue;
			sprintf(fullpath, "%s/%s", path, ent->d_name);
			parse_config(fullpath);
		}
		closedir(dir);
	}

	sprintf(path, "%s/.rocpack/config", home);
	parse_config(path);
}

/* full input file can have all sections */

void parse_input_file(const char* filename)
{
	/* open input file */
	if ((input = fopen(filename, "r")) == NULL)
		errx(EX_NOINPUT, "%s", filename);

	line = 0; shift(); /* load first token */

	parse_settings();

	while(accept(IMPORT)) {
		if (accept('<')) {
			char *filename = get_identifier();
			expect('>');

			import(filename);

			free(filename);
		} else if (ctoken.type == STRING) {
			import(ctoken.value.str);
		} else {
			syntax_error("cannot import file");
		}
	}

	parse_settings();

	while (accept(SHAPE))
		parse_new_shape();

	parse_settings();

	parse_boundary();

	header = line;

	if (ctoken.type == CREATE)
		parse_init_bare();
	else
		parse_init_pack();

	if (particle.size() == 0) {
		syntax_error("no particles were created");
	}

	fclose(input);
}

void write_output_file()
{
	int c, n = 0;
	FILE* output;

	if (output_fname == NULL)
		return;

	if ((input = fopen(input_fname, "r")) == NULL)
		err(1, "error opening file %s", input_fname);

	if ((output = fopen(output_fname, "w")) == NULL)
		err(1, "error opening file %s", output_fname);

	/* set the output scaling factor */

	if (!settings::isset("scaling")) {
		float scaling = settings::get<float>("max_size") / t_curr;
		fprintf(output, "\nset scaling = %f;\n\n", scaling);
	}

	/* copy header */

	while (n < header) {
		c = fgetc(input);
		fputc(c, output);
		if (c == '\n')
			n++;
	}

	/* write particles */

	for (unsigned int i = 0; i < particle.size(); i++) {
		Point x = particle[i]->position();
		Quaternion q = particle[i]->orientation();
		Vector axis = Vector(q[0],q[1],q[2]).normalized();
		float angle = 2.0 * deg(acos(q[3]));
		float size = t_curr * particle[i]->growth_rate();

		float r, g, b, a; particle[i]->get_color(r, g, b, a);

		fprintf(output, "\n%s {\n", particle[i]->shape()->name());
		fprintf(output, "\ttranslate <%11.8f, %11.8f, %11.8f>\n",x[0], x[1], x[2]);
		fprintf(output, "\trotate <%11.8f, %11.8f, %11.8f> %13.8f\n",
			axis[0], axis[1], axis[2], angle < 180.0 ? angle : angle - 360.0);
		fprintf(output, "\tscale %10.8f color %f %f %f tag %u\n",
			size, r, g, b, particle[i]->tag());
		fprintf(output, "}\n");
	}

	fclose(input);
	fclose(output);

#if 0
	/* write contact network */

	FILE* contacts;

	if ((contacts = fopen("contacts.out", "w")) == NULL)
		err(1, "error opening file contacts.out");

	for (unsigned int i = 0; i < particle.size(); i++) {
		contact_network(*particle[i], contacts);
	}

	fclose(contacts);
#endif
}

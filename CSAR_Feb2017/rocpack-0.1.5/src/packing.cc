#include <err.h>
#include <time.h>
#include <cstdio>
#include <vector>

#include "stime.h"
#include "queue.h"
#include "particle.h"
#include "boundary.h"
#include "settings.h"
#include "timer.h"
#include "hgrid.h"
#include "gjk.h"

#include "opengl.h"

extern float target_fraction;
extern int verbose;
#ifdef HAVE_OPENGL
extern int draw, draw_final;
#endif
extern char* output_fname;

static Timer<milliseconds> main_timer;
static double msec = 0, msec_prev = 0, dmsec = 0;

double msec_max;
float rate_min;

void sync()
{
	for (unsigned int i = 0; i < particle.size(); i++) {
		particle[i]->sync(t_curr);
	}
}

//#define TIME_BENCH

inline void progress()
{
	static float f_prev = 1.0, f = 0.0;
	static unsigned long int collisions = 0;

	collisions++; f = pow(t_curr/t_stop, 3.0);

	msec  = main_timer.elapsed();
	dmsec = msec - msec_prev;

	if (f - f_prev > 0.0001 || dmsec > 999 || t_curr >= t_stop) {
		float packing_rate = target_fraction * (f - f_prev) / (dmsec/60000.0);

		if (packing_rate > 9.0) packing_rate = 9.999;

#ifdef EDMC
		if (verbose) {
			fprintf(stderr, "\r%5.2f%% %.4f %.4f/min %.1e ev/s %.1f s   ",
				100*f, target_fraction*f, packing_rate, collisions/(dmsec/1000.0), msec/1000.0);
#ifdef TIME_BENCH
			fprintf(stdout, "%10.4f %.6f\n", msec/1000.0, target_fraction*f);
#endif
		}

#else
		float collision_frequency = collisions / (t_curr - t_prev);

		if (verbose)
			fprintf(stderr, "\r%5.2f%% %.6f %.6f/min %.1e coll/s %5.2f kev/s %.1f s      ",
				100*f, target_fraction*f, packing_rate,
				collision_frequency, (float)(collisions)/dmsec, msec/1000.0);

		if (packing_rate > 0.0 && packing_rate < rate_min) {
			target_fraction *= pow(t_curr/t_stop, 3.0);
			t_stop = t_curr;
		}
#endif

		f_prev = f; collisions = 0;
	}

	if (msec >= msec_max) {
		target_fraction *= pow(t_curr/t_stop, 3.0);
		t_stop = t_curr;
	}
}

void post_run()
{
#ifdef HAVE_OPENGL
	if (draw) {
		return;
	}
#endif

	sync();

	if (verbose)
		progress();

#ifdef HAVE_OPENGL
	if (draw_final) {
		draw_final = false;
		opengl_init(0, NULL);
		glutIdleFunc(render);
		glutMainLoop();
	} else {
		if (verbose)
			fprintf(stderr, "\n");
	}
#else
	if (verbose)
		fprintf(stderr, "\n");
#endif
}

void new_particle_velocities()
{
	for (unsigned int i = 0; i < particle.size(); i++) {
		particle[i]->set_velocity(1.0*Vector::random());
		particle[i]->set_angular_velocity(1.0*Vector::random());
	}
}

void move_back(float dt)
{
	float time = t_curr - dt;

	for (unsigned int i = 0; i < particle.size(); i++) {
		particle[i]->set_time(time);
	}

	t_curr = time; t_prev = 0.0;
}

void stop_particles()
{
	for (unsigned int i = 0; i < particle.size(); i++) {
		particle[i]->set_velocity(Vector(0.0, 0.0, 0.0));
		particle[i]->set_angular_velocity(Vector(0.0, 0.0, 0.0));
	}
}

#ifdef HAVE_OPENGL
void proceed_opengl(){
#ifdef DEBUG_PROGRESS
	/* add blank line to distinguish events */
	fprintf(stderr, "\n"); progress(); fprintf(stderr, "\n");
#endif

	t_prev = t_curr;

	if (t_next > t_stop) {
		t_curr = t_stop; t_next = t_stop;
		sync(); progress(); post_run();
	}

#ifdef SMOOTH_RENDERING
	static const float dt = 0.001;
	while ((t_curr + dt) < t_next - EPSILON) {
		t_curr += dt; sync(); render();
	}
#endif

	t_curr = t_next; sync(); render();

#ifdef DEBUG_EVENTS
	event_queue->next_event()->print();
#else
	static unsigned int n = 0;

	if (n % particle.size()){
		n = 1; progress();
	} else n++;
#endif

	event_queue->next_event()->process();
	event_queue->update();
	t_next = event_queue->next_event()->time();
}
#endif

inline void proceed(){
#ifdef DEBUG_PROGRESS
	/* add blank line to distinguish events */
	fprintf(stderr, "\n"); progress(); fprintf(stderr, "\n");
#endif

	t_prev = t_curr;
	t_curr = t_next < t_stop ? t_next : t_stop;

	if (t_curr >= t_stop) return;

#ifdef DEBUG_EVENTS
	event_queue->next_event()->print();
#else
	static unsigned int n = 0;

	if (n % particle.size()){
		n = 1; progress();
	} else n++;
#endif

	event_queue->next_event()->process();
	event_queue->update();
	t_next = event_queue->next_event()->time();
}

void init()
{
	/* pack initialization */

	float volume = 0.0;
	float min_size = settings::get<float>("min_size");
	float max_size = settings::get<float>("max_size");
	float growth_rate = settings::get<float>("growth_rate");

	assert(min_size > 0.0);
	assert(max_size >= min_size);
	assert(growth_rate > 0.1 && growth_rate < 10.0);
	assert(boundary != NULL);

	if (settings::get<bool>("noinit")) {
		for (unsigned int i = 0; i < particle.size(); i++) {
			float r = growth_rate * particle[i]->growth_rate() / max_size;
			particle[i]->set_growth_rate(r);
			volume += particle[i]->volume();
		}
	} else {
		for (unsigned int i = 0; i < particle.size(); i++) {
			float r = growth_rate * particle[i]->growth_rate() / max_size;
			particle[i]->set_growth_rate(r);
			particle[i]->set_position(boundary->get_random_position());
#ifndef LATTICE_PACKING
			particle[i]->set_orientation(Quaternion::random());
#endif
			volume += particle[i]->volume();
		}
	}

	/* read packing fraction from input file */
	if (target_fraction == 1.0 && settings::isset("packing_fraction"))
		target_fraction = settings::get<float>("packing_fraction");

	t_start = t_prev = t_curr = 0.0;
	t_stop  = pow(target_fraction * boundary->volume() / volume, 1.0/3.0);

	int levels = 0;

	if (settings::isset("hgrid_levels")) {
		levels = settings::get<int>("hgrid_levels");
	} else {
		while (max_size/min_size > (2 << levels)) levels++;

		if (settings::isset("max_hgrid_levels")) {
			set_max(levels, settings::get<int>("max_hgrid_levels"));
		}
	}

	hgrid = new HGrid(4.01 * growth_rate * t_stop, levels);

	for (unsigned int i = 0; i < particle.size(); i++) {
		particle[i]->set_event_id(i);
		particle[i]->set_hash(0xffffffffffffff); /* necessary! */
		hgrid->insert(particle[i]);
	}

	event_queue = new EventQueue(particle.size());
}

#ifdef SPHERE_PRE_PACK

/* pre-pack particles using bounding spheres */

void main_loop_spheres()
{
	float t_old = 0.0;
	unsigned int n = 0;
	unsigned int j = 1;

	std::vector<Shape*> shapes;
	Shape* sphere = create_shape("sphere");

	/* temporarily change shapes to spheres */

	for (unsigned int i = 0; i < particle.size(); i++) {
		shapes.push_back(particle[i]->shape());
		particle[i]->set_shape(sphere);
	}

	/* prepare event queue */

	for (unsigned int i = 0; i < particle.size(); i++) {
		if (verbose)
			fprintf(stderr, "\rinit [%3lu%%]", 100 * (i+1) / particle.size());
		Event *event = new Event(i, particle[i]);
		event_queue->push(event);
	}

	t_next = event_queue->next_event()->time();

	for(;;) {
		t_prev = t_curr;
		t_curr = t_next < t_stop ? t_next : t_stop;

		if (t_curr >= t_stop) return;

		if (verbose) {
			if (n % particle.size()){
				n = 1; progress();
			} else n++;
		}

		if (j++ > 10000) {
			j = 1;
			if ((t_curr - t_old)/t_stop < 0.001)
				break;
			else
				t_old = t_curr;
		}

		event_queue->next_event()->process();

		event_queue->update();

		t_next = event_queue->next_event()->time();
	}

	/* restore original shapes */

	for (unsigned int i = 0; i < particle.size(); i++)
		particle[i]->set_shape(shapes[i]);

	sync();

	stop_particles();

	t_curr = 0.0;

	/* recompute events with original shapes */
	event_queue->update_all();
}

#endif

void main_loop()
{
	main_timer.start();

	if (settings::isset("max_runtime"))
		msec_max = 1000 * settings::get<double>("max_runtime");
	else
		msec_max = DBL_MAX;

	if (settings::isset("packing_rate_min"))
		rate_min = settings::get<float>("packing_rate_min");
	else
		rate_min = 0.0; /* never stop due to packing rate */

	init();

#ifdef SPHERE_PRE_PACK
	main_loop_spheres();
#else
	for (unsigned int i = 0; i < particle.size(); i++) {
		if (verbose)
			fprintf(stderr, "\rinit [%3lu%%]", 100 * (i+1) / particle.size());
		Event *event = new Event(i, particle[i]);
		event_queue->push(event);
	}
#endif

	t_next = event_queue->next_event()->time();

#ifdef HAVE_OPENGL
	if (draw) {
		opengl_init(0, NULL);
		glutIdleFunc(proceed_opengl);
		glutMainLoop();
	} else {
		while(t_curr < t_stop)
			proceed();

		post_run();
	}
#else
	while(t_curr < t_stop)
		proceed();

	post_run();
#endif

	if (!verbose) {
		msec = main_timer.elapsed();

		fprintf(stdout, "%.4f %.1f\n", target_fraction * pow(t_curr/t_stop, 3.0), msec/1000.0);
	}
}


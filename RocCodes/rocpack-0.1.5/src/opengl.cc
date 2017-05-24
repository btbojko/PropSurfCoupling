#include "opengl.h"

#ifdef HAVE_OPENGL
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef HAVE_GD
#include "gd.h"
#endif

#include "hgrid.h"
#include "queue.h"
#include "particle.h"
#include "boundary.h"
#include "settings.h"
#include "collision.h"
#include "stime.h"
#include "convex.h"
#include "gjk.h"

extern float t_start, t_curr, t_stop;

GLfloat light_color[][4] = {
	{ 0.5f, 0.5f, 0.5f, 1.0f},
	{ 0.5f, 0.5f, 0.5f, 1.0f}
};

GLfloat light_position[][4] = {
	{ 5.0f, 5.0f, 5.0f, 1.0f},
	{-5.0f, 5.0f, 5.0f, 1.0f}
};

GLfloat ambient_color[] = {0.5f, 0.5f, 0.5f, 1.0f};

static int window;
static int width  = 1024;
static int height = 1024;

static float zoom = 2.1;
static bool in_motion = false;
int last_x = 0, last_y = 0, curr_x = 0, curr_y = 0;
static Quaternion orientation(0.228533, 0.242249, 0.036453, 0.942209); // 0.0, 0.0, 0.0, 1.0);

int save_screenshot(const char*);

void opengl_keyboard(unsigned char key, int x, int y) {
	switch (key) {
		case 's':
		case 'p':
			save_screenshot("screenshot.png");
			break;

		case 'i':
		case 'z':
		case '+':
			zoom -= 0.05; break;

		case 'o':
		case 'Z':
		case '-':
			zoom += 0.05; break;

		case  27: /* Escape */
		case 'q':
		case 'Q':
			fprintf(stdout,"\n");
#ifdef __APPLE__
			/* Apple's GLUT hangs if we just destroy the window... */
			exit(0);
#else
			glutDestroyWindow(window);
#endif
	}
}

void opengl_mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		in_motion = true;
		last_x = curr_x = x;
		last_y = curr_y = y;
	} else {
		in_motion = false;
	}
}

void opengl_motion(int x, int y)
{
	if (in_motion) { curr_x = x; curr_y = y; }
}

Vector trackball_vector(int x, int y)
{
	Vector v(2.0*x/width - 1.0, -(2.0*y/height - 1.0), 0.0);

	if (norm2(v) <= 0.5)
		v[2] = sqrt(1.0 - norm2(v));
	else
		v[2] = 0.5 / norm(v);

	return v.normalized();
}

void compute_camera_rotation(Vector& axis, float& angle)
{
	Vector v1 = trackball_vector(last_x, last_y);
	Vector v2 = trackball_vector(curr_x, curr_y);

	axis  = (v1 ^ v2).normalized();
	angle = acos(v1*v2);
}

void opengl_reshape_window(int w, int h){
	width = w; height = h;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (float) width / (float) height, 0.1, 50.0);
}

void opengl_init(int argc, char* argv[]){

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	if (settings::isset("window_size")) {
		width  = settings::get<int>("window_size");
		height = settings::get<int>("window_size");
	}

	if (settings::isset("window_width"))
		width  = settings::get<int>("window_width");

	if (settings::isset("window_height"))
		height = settings::get<int>("window_height");

	glutInitWindowSize(width, height);

	window = glutCreateWindow("Rocpack");

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
    glEnable(GL_NORMALIZE); /* important to make shading work */
    //glEnable(GL_BLEND);     /* important to make transparency work */
	//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel(GL_SMOOTH);

	glutKeyboardFunc(opengl_keyboard);
	glutMouseFunc(opengl_mouse);
	glutMotionFunc(opengl_motion);
	glutReshapeFunc(opengl_reshape_window);
	glutDisplayFunc(render);
}

void render(){
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

#if 0
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_color);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color[0]);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position[0]);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_color[1]);
	glLightfv(GL_LIGHT1, GL_POSITION, light_position[1]);
#endif

	glTranslatef(0.0, 0.0, -zoom);

	if (last_x != curr_x || last_y != curr_y) {
		Vector axis; float angle;
		compute_camera_rotation(axis, angle);
		axis = orientation.conjugated() ^ axis ^ orientation;
		orientation ^= Quaternion(axis, angle);
		angle = 2.0 * deg(acos(orientation[3]));
		glRotatef(angle, orientation[0], orientation[1], orientation[2]);
		last_x = curr_x; last_y = curr_y;
	} else {
		float angle = 2.0 * deg(acos(orientation[3]));
		glRotatef(angle, orientation[0], orientation[1], orientation[2]);
	}

	glColor3f(0.0000, 0.0000, 0.0000);

	boundary->draw();

#ifdef DEBUG_DRAW
	hgrid->draw();
#endif

	/* hilight particles in debug mode */

#ifndef DEBUG_DRAW
	for(unsigned int i = 0; i < particle.size(); i++){ particle[i]->draw(); }
#else
	int id1 = event_queue->next_event()->id();
	int id2 = event_queue->next_event()->secondary_id();

	if (id2 == -1) { /* only one particle in event */

		/* hilight current event's particle */

		glColor3f(0.7843, 0.4941, 0.1529);
		particle[id1]->draw();

		/* draw remaining particles */

		glColor3f(0.1529, 0.4941, 0.7843);
		for(unsigned int i = 0; i < (unsigned int)(id1); i++){ particle[i]->draw(); }
		for(unsigned int i = id1+1; i < particle.size(); i++){ particle[i]->draw(); }

	} else {
		unsigned int hl0 = min<unsigned int>(id1, id2);
		unsigned int hl1 = max<unsigned int>(id1, id2);

		/* hilight colliding particles */

		glColor3f(0.7843, 0.4941, 0.1529);

		particle[hl0]->draw();
		particle[hl1]->draw();

		/* draw remaining particles */

		glColor3f(0.1529, 0.4941, 0.7843);
		for(unsigned int i = 0; i < hl0; i++){ particle[i]->draw(); }
		for(unsigned int i = hl0+1; i < hl1; i++){ particle[i]->draw(); }
		for(unsigned int i = hl1+1; i < particle.size(); i++){ particle[i]->draw(); }

#ifdef DEBUG_CLOSEST_POINTS
		/* draw closest points of colliding particles */

		Point pA, pB;

		if (intersect(*particle[hl0], *particle[hl1], t_curr)) {
			fprintf(stderr, "error: intersection found at t = %.20f: particle %d with %d\n", t_curr, hl0, hl1);
			fprintf(stderr, "distance: %.20f\n", distance(*particle[hl0], *particle[hl1], t_curr));
			exit(1);
		}

		closest_points(*particle[hl0], *particle[hl1], pA, pB, t_curr);

		glColor3f(0.0000, 0.0000, 0.0000);

		glBegin(GL_LINES);
		glVertex3f(pA[0], pA[1], pA[2]);
		glVertex3f(pB[0], pB[1], pB[2]);
		glEnd();
#endif
	}
#endif

	glutSwapBuffers();
}

#ifdef HAVE_GD
int save_screenshot(const char *filename){
    FILE *png;
    GLubyte *OpenGLimage, *p;
    gdImagePtr image;
    unsigned int r, g, b;
    int i, j, rgb;

    png = fopen(filename, "wb");

    if (png == NULL) {
        fprintf(stderr, "warning: unable to write to %s\n", filename);
        return 1;
    }

    OpenGLimage = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

    if(OpenGLimage == NULL){
        fprintf(stderr, "error allocating image:%s\n", filename);
        return 2;
    }

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, OpenGLimage);
    p = OpenGLimage;
    image = gdImageCreateTrueColor(width,height);

    for (i = height-1 ; i>=0; i--) {
        for(j=0;j<width;j++){
                r=*p++; g=*p++; b=*p++;
                rgb = (r<<16)|(g<<8)|b;
                gdImageSetPixel(image,j,i,rgb);
        }
    }

    gdImagePng(image,png);
    fclose(png);
    gdImageDestroy(image);

	return 0;
}
#else
int save_screenshot(const char *filename) { return 0; }
#endif

#endif

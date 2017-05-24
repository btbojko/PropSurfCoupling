#ifndef _OPENGL_H_
#define _OPENGL_H_

#include "config.h"

#ifdef HAVE_OPENGL

#if HAVE_WINDOWS_H && defined(_WIN32)
# include <windows.h>
#endif
#if defined(HAVE_GL_GLUT_H)
# include <GL/glut.h>
#elif defined(HAVE_GLUT_GLUT_H)
# include <GLUT/glut.h>
#else
# error no glut.h
#endif

void opengl_init(int argc, char* argv[]);

void render();

#endif

#endif

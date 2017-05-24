#ifndef _DEBUG_H
#define _DEBUG_H

#include <cassert>

#ifdef NDEBUG

#define DEBUG( ... ) /* nothing */
#define DEBUG_INT(n) /* nothing */
#define DEBUG_FLT(x) /* nothing */
#define DEBUG_VEC(v) /* nothing */

#else

#include <cstdio>
#include <cstdarg>

#ifdef DEBUG
#undef DEBUG
#endif

#define DEBUG_LOC fprintf(stderr,"%s:%d ", __FILE__, __LINE__)
#define DEBUG( ... ) do { DEBUG_LOC; fprintf(stderr, __VA_ARGS__ ); } while (0)
#define DEBUG_INT(n) fprintf(stderr, "%s = %d\n", #n, n)
#define DEBUG_FLT(x) fprintf(stderr, "%s = %f\n", #x, x)
#define DEBUG_VEC(v) fprintf(stderr, "%s = [%f; %f %f] (%f)\n", \
                               #v, v[0], v[1], v[2], norm(v))
#endif

#endif

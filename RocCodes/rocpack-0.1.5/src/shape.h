#ifndef _SHAPE_H
#define _SHAPE_H

#include "scalar.h"

class Point;
class Vector;
class Matrix;
class Transform;

enum ShapeType  { SPHERE, CONVEX, CONCAVE, COMPOUND };

class Shape {
  public:
	virtual ~Shape() {}
	virtual ShapeType type() const = 0;
	virtual const char* name() const = 0;
	virtual float bounding_radius() const = 0;
	virtual Matrix inertia() const = 0;
	virtual Matrix inverse_inertia() const = 0;
	virtual float volume() const = 0;
#ifdef HAVE_OPENGL
	virtual void draw() const = 0;
#endif
};

Shape* create_shape(const char* shape);
bool register_shape(const char *shape, Shape* shape_ptr);

#endif

#ifndef _CONVEX_H
#define _CONVEX_H

#include "shape.h"

class Convex : public Shape {
  public:
	virtual ~Convex() {}
	virtual ShapeType type() const { return CONVEX; }
	virtual Point support(const Vector& v) const = 0;
};

#endif

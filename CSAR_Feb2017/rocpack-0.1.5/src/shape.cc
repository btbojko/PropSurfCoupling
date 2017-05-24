#include <map>
#include <string>

#include <err.h>
#include <sysexits.h>

#include "shape.h"

static std::map<std::string, Shape*> shape_map;

bool register_shape(const char *shape, Shape* shape_ptr)
{
	return (shape_map[shape] = shape_ptr) != NULL;
}

Shape* create_shape(const char* shape)
{
	if (shape_map.find(shape) == shape_map.end())
		errx(EX_SOFTWARE, "unknown shape: %s", shape);

	return shape_map[shape];
}

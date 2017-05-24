ACLOCAL_AMFLAGS = -I autotools/m4

install-exec-home:
	make bindir=$(HOME)/bin install-exec

install-shapes-home:
	make shapesdir=$(HOME)/.$(PACKAGE)/shapes install-dist_shapesDATA

install-home: install-exec-home install-shapes-home

AM_CXXFLAGS = -O3 -march=native

bin_PROGRAMS = pack pack-ls rdf pack2mesh pack2bin

rdf_SOURCES = src/rdf.cc
rdf_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

pack_CPPFLAGS = -DEDMC
pack_ls_CPPFLAGS = -UEDMC

pack_LDADD = $(OPENGL_LIBS) $(GD_LIBS)
pack_ls_LDADD = $(OPENGL_LIBS) $(GD_LIBS)

pack2mesh_LDADD = $(OPENGL_LIBS) $(GD_LIBS) $(TIFF_LIBS)
pack2bin_LDADD = $(OPENGL_LIBS) $(GD_LIBS) $(TIFF_LIBS)

dist_bin_SCRIPTS = \
	scripts/obj2shape \
	scripts/off2shape \
	scripts/pack2rdf  \
	scripts/shapegen

dist_doc_DATA = doc/manual.pdf doc/rocpack.pdf

shapesdir = $(pkgdatadir)/shapes

dist_shapes_DATA = \
	shapes/adn \
	shapes/antiprism4 \
	shapes/antiprism5 \
	shapes/antiprism6 \
	shapes/cube \
	shapes/cuboctahedron \
	shapes/dodecahedron \
	shapes/ellipsoids20x20 \
	shapes/hmx \
	shapes/icosahedron \
	shapes/icosidodecahedron \
	shapes/octahedron \
	shapes/petn \
	shapes/prism3 \
	shapes/prism5 \
	shapes/prism6 \
	shapes/pyramid4 \
	shapes/pyramid5 \
	shapes/pyramid6 \
	shapes/reflected_snub_cube \
	shapes/reflected_snub_dodecahedron \
	shapes/rhombicosidodecahedron \
	shapes/rhombicuboctahedron \
	shapes/rounded_cube \
	shapes/snub_cube \
	shapes/snub_dodecahedron \
	shapes/snub_prism3 \
	shapes/square_rod \
	shapes/tetrahedron \
	shapes/truncated_cube \
	shapes/truncated_cuboctahedron \
	shapes/truncated_dodecahedron \
	shapes/truncated_icosahedron \
	shapes/truncated_icosidodecahedron \
	shapes/truncated_octahedron \
	shapes/truncated_tetrahedron

if EXAMPLES
examplesdir = $(pkgdatadir)/examples

dist_examples_SCRIPTS = \
	examples/ap-propellant \
	examples/boundary-types \
	examples/crystals \
	examples/cylinders \
	examples/disks \
	examples/ellipsoids \
	examples/icosahedra \
	examples/new-shape \
	examples/polydisperse \
	examples/simple \
	examples/speed-test \
	examples/visual

dist_examples_DATA = examples/readbin.f90
endif

pack_SOURCES = \
	src/boundary.h \
	src/collision.h \
	src/convex.h \
	src/cylinder.h \
	src/debug.h \
	src/ellipsoid.h \
	src/event.h \
	src/gjk.h \
	src/hgrid.h \
	src/matrix.h \
	src/opengl.h \
	src/packing.h \
	src/parse.h \
	src/particle.h \
	src/point.h \
	src/polyhedron.h \
	src/quaternion.h \
	src/queue.h \
	src/random.h \
	src/scalar.h \
	src/settings.h \
	src/shape.h \
	src/simplex.h \
	src/sphere.h \
	src/stack.h \
	src/state.h \
	src/stime.h \
	src/timer.h \
	src/token.h \
	src/transform.h \
	src/vector.h \
	src/vmath.h \
	src/boundary.cc \
	src/collision.cc \
	src/cylinder.cc \
	src/ellipsoid.cc \
	src/event.cc \
	src/gjk.cc \
	src/hgrid.cc \
	src/main.cc \
	src/opengl.cc \
	src/packing.cc \
	src/parse.cc \
	src/particle.cc \
	src/polyhedron.cc \
	src/queue.cc \
	src/settings.cc \
	src/shape.cc \
	src/simplex.cc \
	src/sphere.cc \
	src/token.cc

pack_ls_SOURCES = \
	src/boundary.h \
	src/collision.h \
	src/convex.h \
	src/cylinder.h \
	src/debug.h \
	src/ellipsoid.h \
	src/event.h \
	src/gjk.h \
	src/hgrid.h \
	src/matrix.h \
	src/opengl.h \
	src/packing.h \
	src/parse.h \
	src/particle.h \
	src/point.h \
	src/polyhedron.h \
	src/quaternion.h \
	src/queue.h \
	src/random.h \
	src/scalar.h \
	src/settings.h \
	src/shape.h \
	src/simplex.h \
	src/sphere.h \
	src/stack.h \
	src/state.h \
	src/stime.h \
	src/timer.h \
	src/token.h \
	src/transform.h \
	src/vector.h \
	src/vmath.h \
	src/boundary.cc \
	src/collision.cc \
	src/cylinder.cc \
	src/ellipsoid.cc \
	src/event.cc \
	src/gjk.cc \
	src/hgrid.cc \
	src/main.cc \
	src/opengl.cc \
	src/packing.cc \
	src/parse.cc \
	src/particle.cc \
	src/polyhedron.cc \
	src/queue.cc \
	src/settings.cc \
	src/shape.cc \
	src/simplex.cc \
	src/sphere.cc \
	src/token.cc

pack2mesh_SOURCES = \
	src/boundary.h \
	src/collision.h \
	src/convex.h \
	src/cylinder.h \
	src/debug.h \
	src/ellipsoid.h \
	src/event.h \
	src/gjk.h \
	src/hgrid.h \
	src/matrix.h \
	src/opengl.h \
	src/packing.h \
	src/parse.h \
	src/particle.h \
	src/point.h \
	src/polyhedron.h \
	src/quaternion.h \
	src/queue.h \
	src/random.h \
	src/scalar.h \
	src/settings.h \
	src/shape.h \
	src/simplex.h \
	src/sphere.h \
	src/stack.h \
	src/state.h \
	src/stime.h \
	src/timer.h \
	src/token.h \
	src/transform.h \
	src/vector.h \
	src/vmath.h \
	src/boundary.cc \
	src/collision.cc \
	src/cylinder.cc \
	src/ellipsoid.cc \
	src/event.cc \
	src/gjk.cc \
	src/hgrid.cc \
	src/opengl.cc \
	src/pack2mesh.cc \
	src/packing.cc \
	src/parse.cc \
	src/particle.cc \
	src/polyhedron.cc \
	src/queue.cc \
	src/settings.cc \
	src/shape.cc \
	src/simplex.cc \
	src/sphere.cc \
	src/token.cc

pack2bin_SOURCES = \
	src/boundary.h \
	src/collision.h \
	src/convex.h \
	src/cylinder.h \
	src/debug.h \
	src/ellipsoid.h \
	src/event.h \
	src/gjk.h \
	src/hgrid.h \
	src/matrix.h \
	src/opengl.h \
	src/packing.h \
	src/parse.h \
	src/particle.h \
	src/point.h \
	src/polyhedron.h \
	src/quaternion.h \
	src/queue.h \
	src/random.h \
	src/scalar.h \
	src/settings.h \
	src/shape.h \
	src/simplex.h \
	src/sphere.h \
	src/stack.h \
	src/state.h \
	src/stime.h \
	src/timer.h \
	src/token.h \
	src/transform.h \
	src/vector.h \
	src/vmath.h \
	src/boundary.cc \
	src/collision.cc \
	src/cylinder.cc \
	src/ellipsoid.cc \
	src/event.cc \
	src/gjk.cc \
	src/hgrid.cc \
	src/opengl.cc \
	src/pack2bin.cc \
	src/packing.cc \
	src/parse.cc \
	src/particle.cc \
	src/polyhedron.cc \
	src/queue.cc \
	src/settings.cc \
	src/shape.cc \
	src/simplex.cc \
	src/sphere.cc \
	src/token.cc

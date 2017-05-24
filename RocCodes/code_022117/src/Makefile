SHELL		= /bin/sh
.SUFFIXES:	.o .f90

OSTYPE		= $(shell uname)

EXEC 		= rocfire3d

OBJS            = data_types.o data_str.o mg_solver.o lusolver.o restart.o \
		 driver.o grid.o ic.o\
                pressure.o rate.o bcnew.o velocity.o update_phi.o\
                update_psinew.o pyrolysis.o flux.o conductivity.o\
                txy_2_1.o get_proc_decomp.o get_mg_decomp.o quadrants.o ghost.o global_data.o\
                bc_data.o mpi.o allocation.o qxy_2_1.o\
                bcq.o projection.o\
                solid_prop.o read_pack.o\
                test.o implicit_bc.o explicit_residual.o wgauss_relax.o\
                globalparameter.o viz1.o viz2.o \
		newroutines.o read_matrix.o statistics.o radiation.o



MX		= mod
CD		= \cd
LN		= ln -sf
RM		= rm -f
ECHO            = echo

#----------------- options for Blue Horizon (San Diego)
CXX		= g++
F90C		= mpif90
F77C		= mpif77
F90FLAGS	= -O3 #-Q #-qsuffix=f=f90  -qinitauto -qthreaded -qstrict -qtune=auto
F77FLAGS	= -O3 #-Q #-qinitauto -qthreaded -qstrict -qtune=auto
LD		= $(F90C)
LDLIBS = -lmpich


all :	note  run

note:  
	$(ECHO) "Building $(EXEC) on a $(OSTYPE) system. $(F90C)"
	-$(RM) $(EXEC)

$(EXEC): 	$(OBJS)
	$(LD)  $(LDFLAGS) $(FFLAGS) -o $(EXEC) $(OBJS) $(LDLIBS)

run : 		$(EXEC)

seq : 
	-$(CD) orig_code;$(MAKE)

#-------------------------------- dependencies
velocity.o: velocity.f90 global_data.o mpi.o
restart.o : restart.f90 global_data.o mpi.o statistics.o wgauss_relax.o
data_types.o: data_types.f90
data_str.o: data_str.f90 global_data.o mpi.o data_types.o
mg_solver.o: mg_solver.f90 global_data.o mpi.o data_types.o bc_data.o viz1.o
lusolver.o: lusolver.f90 global_data.o mpi.o
driver.o: driver.f90 global_data.o mpi.o wgauss_relax.o restart.o statistics.o radiation.o
grid.o: grid.f90 global_data.o mpi.o data_types.o
ic.o: ic.f90 global_data.o mpi.o
pressure.o: pressure.f90 global_data.o mpi.o
rate.o: rate.f90 global_data.o mpi.o
bcnew.o: bcnew.f90 global_data.o mpi.o wgauss_relax.o
update_phi.o: update_phi.f90 global_data.o mpi.o restart.o
update_psinew.o:  global_data.o mpi.o
pyrolysis.o:  global_data.o mpi.o
flux.o:  global_data.o mpi.o
conductivity.o:  global_data.o mpi.o
txy_2_1.o:  global_data.o mpi.o wgauss_relax.o
get_proc_decomp.o:  global_data.o mpi.o
get_mg_decomp.o:  global_data.o mpi.o data_str.o
quadrants.o:  global_data.o mpi.o
ghost.o:  global_data.o mpi.o data_types.o
global_data.o:   data_types.o
bc_data.o: global_data.o mpi.o
mpi.o:  global_data.o
allocation.o:  global_data.o mpi.o wgauss_relax.o
qxy_2_1.o:  global_data.o mpi.o
bcq.o:  global_data.o mpi.o
projection.o:  global_data.o mpi.o data_types.o mg_solver.o
solid_prop.o:  global_data.o mpi.o
read_pack.o:  global_data.o mpi.o radiation.o
test.o:  global_data.o mpi.o
implicit_bc.o:  global_data.o mpi.o wgauss_relax.o
explicit_residual.o:  global_data.o mpi.o
wgauss_relax.o:  global_data.o mpi.o data_types.o bc_data.o lusolver.o
globalparameter.o:  global_data.o mpi.o
viz1.o:  global_data.o mpi.o data_types.o statistics.o
viz2.o:  global_data.o mpi.o
newroutines.o:  global_data.o mpi.o data_types.o
read_matrix.o:  global_data.o mpi.o
statistics.o:  global_data.o mpi.o
radiation.o: global_data.o mpi.o  mg_solver.o
#-------------------------------- implicit rules
%.o:  %.f90 
	$(F90C) -c $(F90FLAGS) $< 
%.o:  %.f 
	$(F77C) -c $(F77FLAGS) $< 
.C.o : 
	$(CXX) -c $(OFLAGS)  $<
###%.o: %.f90
###	$(F90C) -c $(F90FLAGS) $<  -o $@
###%.o: %.f
###	$(F77) -c $(FFLAGS) $<  -o $@
#-------------------------------- clean

new :	clean note run 



clean:	runclean
	-$(RM) -f *~ *.o *.mod *.s $(EXEC)
runclean:
	-$(RM) -f core fort.* $(EXEC)

reset:
	rm -f *.dat D/*
#--------------------------------

SHELL		= /bin/sh
.SUFFIXES:	.o .f90 .for

OSTYPE		= $(shell uname)

EXEC 		= rocfire2d
OBJS		=  allocation.o \
		avg_space.o \
 		bc.o \
		bc_data.o \
		bcq.o \
		binder_rate.o \
		conductivity.o \
		data_str.o \
		driver.o \
		explicit_residual.o \
		ed_input.o\
		edwards.o\
		flux.o \
		gauss_relax.o \
		get_mg_decomp.o \
		get_proc_decomp.o \
		ghost.o \
		global_data.o \
		globalparameter.o \
		grid.o \
		ic.o \
		implicit_bc.o \
		lusolver.o  \
		mg_solver.o \
		mpi.o \
		newroutines.o \
		precmod.o\
		pressure.o \
		projection.o \
		psi.o \
		pyrolysis.o \
		rate.o \
		read_pack.o \
		ree2d.o \
		restart.o \
		sandwich_phi.o \
		solid_prop.o  \
		stringmod.o\
		test.o		 \
		txy_2_1.o \
		update_phi.o \
		update_psi.o \
		uvw_2_1.o \
		velocity.o \
		visualization.o \
		wgauss_relax.o \

#echo "Building $(EXEC) for a $(OSTYPE) system."
#--------------------------------

MX		= mod
CD		= \cd
LN		= ln -sf
RM		= rm -f
ECHO            = echo

#----------------- options for Blue Horizon (San Diego)
CXX		= mpicc
F90C		= mpif90
F77C		= mpif77
F90FLAGS   = -O3
#F90FLAGS	= -O3 -Q -pthread
#-O3 -Q -qsuffix=f=f90  -qinitauto -qthreaded -qstrict -qtune=auto
#F77FLAGS	= -O3 -Q -qinitauto -qthreaded -qstrict -qtune=auto
LD		= $(F90C)
#LDFLAGS = -bind_at_load
#LDLIBS = -lpmpich -lmpich -lmpe 
#LDLIBS = -lmpich -lmpe

#-------------------------------- default action
all :	note run

note:  
	$(ECHO) "Building $(EXEC) on a $(OSTYPE) system."
 
$(EXEC): 	$(OBJS) 
	$(LD) -o $(EXEC) $(LDFLAGS) $(OBJS) $(LDLIBS)
run : 		$(EXEC)
seq : 
	-$(CD) orig_code;$(MAKE)

#-------------------------------- dependencies
allocation.o: global_data.o bc_data.o mpi.o
avg_space.o: global_data.o  mpi.o
bc.o:  global_data.o bc_data.o mpi.o
bcq.o:  global_data.o  mpi.o
binder_rate.o : global_data.o
conductivity.o : global_data.o  mpi.o
data_str.o: data_types.o global_data.o  mpi.o
driver.o:  global_data.o bc_data.o mpi.o
explicit_residual.o: global_data.o  mpi.o
flux.o:  global_data.o  mpi.o
gauss_relax.o:  data_types.o global_data.o  mpi.o mg_solver.o 
get_mg_decomp.o:  data_types.o global_data.o  mpi.o 
ghost.o:  data_types.o global_data.o  mpi.o
global_data.o: data_types.o 
globalparameter.o:  global_data.o  mpi.o stringmod.o
grid.o:  data_types.o global_data.o  mpi.o
ic.o: global_data.o  mpi.o
implicit_bc.o: global_data.o  mpi.o
lusolver.o:  global_data.o  mpi.o
mg_solver.o: global_data.o  mpi.o lusolver.o
newroutines.o: global_data.o  mpi.o data_types.o
pressure.o: global_data.o  mpi.o
projection.o: global_data.o  mpi.o data_types.o lusolver.o
psi.o: global_data.o  mpi.o
pyrolysis.o: global_data.o  mpi.o
rate.o: global_data.o  mpi.o
read_pack.o: global_data.o  mpi.o
ree2d.o:  global_data.o  mpi.o
restart.o: global_data.o  mpi.o
sandwich_phi.o: global_data.o  mpi.o
stringmod.o: precmod.o
test.o:  data_types.o global_data.o  mpi.o mg_solver.o
txy_2_1.o:  data_types.o global_data.o  mpi.o mg_solver.o
update_phi.o: global_data.o  mpi.o
uvw_2_1.o:  global_data.o  mpi.o
velocity.o:  global_data.o  mpi.o
visualization.o:  global_data.o  mpi.o
wgauss_relax.o: data_types.o global_data.o  mpi.o bc_data.o lusolver.o

#-------------------------------- implicit rules
%.o: %.f90
	$(F90C) -c $(F90FLAGS) $<
%.o: %.for
	$(F77C) -c $(F77FLAGS) $<
%.o: %.F90
	$(F90B) -c $(F77FLAGS) $<
#-------------------------------- clean
new :	clean note run
clean:
	-$(RM) -f *~ *.o *.mod *.s *.ini *.vtk *.dat ED_INPUT.DAT $(EXEC)
runclean:
	-$(CD) data;$(RM) -f core fort.* $(EXEC)
cleanall: runclean
	-$(CD) orig_code;$(MAKE) clean
#--------------------------------


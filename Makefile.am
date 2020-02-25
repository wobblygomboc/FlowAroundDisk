# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= \
flow_around_disk

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
flow_around_disk_SOURCES = flow_around_disk.cc tetmesh_faceted_surfaces.h

# Required libraries: 
flow_around_disk_LDADD = -L@libdir@  -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@  

#EXTRA_DIST = boundaries_2d.lay boundaries_2d.png convergence_2d.lay run2d.bash




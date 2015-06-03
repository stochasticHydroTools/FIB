######################################################################
## Here specify the location of the IBAMR source, the FIB source, and the location
## where IBAMR has been built.
FIB_SRC_DIR = $(HOME)/FluctHydro/FIB
IBAMR_SRC_DIR = $(HOME)/FluctHydro/git_ibamr

# This should be the directory with the IBAMR object files.
IBAMR_BUILD_DIR = /fluct/delong/ibamr-opt


######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc
VPATH=$(FIB_SRC_DIR)

## Update compiler flags.
LIBS     += -L/usr/lib64/ -lm
CFLAGS   +=  -O3 -Wall 
CXXFLAGS +=  -isystem -Wall 
CPPFLAGS += -Wextra -pthread

SOURCES = main.C IBBrownianBlobHierarchyIntegrator.C IBBrownianBlobHierarchyIntegrator.h NonbondedForceEvaluator.C NonbondedForceEvaluator.h WallForceEvaluator.C WallForceEvaluator.h Wall.C Wall.h
OBJS = main.o IBBrownianBlobHierarchyIntegrator.o NonbondedForceEvaluator.o WallForceEvaluator.o Wall.o

#################################################################
default:
	@echo "make one of: main2d, main3d"

main2d:
	if (test -f stamp-3d); then $(MAKE) clean; fi
	touch stamp-2d
	$(MAKE) PDIM=2 main-2d

main3d:
	if (test -f stamp-2d); then $(MAKE) clean; fi
	touch stamp-3d
	$(MAKE) PDIM=3 main-3d

main-2d: $(IBAMR_LIB_2D) $(IBTK_LIB_2D) $(OBJS) $(SOURCES)
	$(CXX)  $(CXXFLAGS)  $(LDFLAGS) $(OBJS) \
	$(IBAMR_LIB_2D) $(IBTK_LIB_2D) $(LIBS) -DNDIM=$(PDIM) -o main2d

main-3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) $(SOURCES) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) \
	$(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LIBS) -DNDIM=$(PDIM) -o main3d

clean:
	$(RM) main2d main3d
	$(RM) *.o *.lo *.objs *.ii *.int.c stamp-[23]d
	$(RM) -r .libs

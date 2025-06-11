# Makefile for CLM surface data processing

include ./Makeoptions

# name of executable
EXECUTABLE = mkgrd.x

# Source directory
SRC_DIR = ./src

# Build directory
BLD_DIR = .bld

# Set VPATH to search for source files, not used in current version
VPATH = $(SRC_DIR)

####################################################################

.PHONY: all mkdir_build clean

# Default target
all: $(EXECUTABLE)

mkdir_build :
	mkdir -p $(BLD_DIR)

# Define object files with build directory
OBJS = $(addprefix $(BLD_DIR)/, \
	consts_coms.o \
	MOD_utilities.o \
	icosahedron.o \
	blas.o \
	lapack.o \
	MOD_data_preprocessing.o \
	MOD_grid_preprocessing.o \
	MOD_mask_processing.o \
	MOD_Area_judge.o \
	MOD_GetContain.o \
	MOD_GetRef.o \
	MOD_refine.o \
	MOD_mask_postprocessing.o \
	MOD_namelist.o \
    MOD_voronoi_geometry.o \
	MOD_grid_initialization.o \
	mkgrd.o \
)

		
####################################################################

# Create symbolic link to executable in current directory
$(EXECUTABLE): $(BLD_DIR)/${EXECUTABLE}
	ln -sf $(BLD_DIR)/${EXECUTABLE} $@
	@echo 'program for making unstructure mesh surface data have been compiled successfully !'

$(BLD_DIR)/${EXECUTABLE} : mkdir_build $(OBJS)
	${FF} ${FOPTS} $(OBJS) -o $@ ${LDFLAGS}

$(BLD_DIR)/%.o : %.F90 mkdir_build
	${FF} -c ${FOPTS} -module $(BLD_DIR) $(INCLUDE_DIR) -o $@ $< 

clean :
	${RM} -f $(BLD_DIR)/*.o $(BLD_DIR)/*.mod $(BLD_DIR)/${EXECUTABLE} ${EXECUTABLE}

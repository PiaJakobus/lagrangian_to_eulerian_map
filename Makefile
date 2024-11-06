# Compiler and flags
FC = gfortran -fopenmp
flag = -O0 -fbounds-check -mcmodel=large -r16 -g -fbacktrace -fbounds-check \
			 -ffpe-trap=zero,overflow,underflow -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter 


# Library paths (adjust these if the libraries are in non-standard locations)
#LAPACK_PATH = /usr/lib/x86_64-linux-gnu   # Adjust if LAPACK is in a non-standard location
#BLAS_PATH = /usr/lib/x86_64-linux-gnu     # Adjust if BLAS is in a non-standard location

# LAPACK and BLAS libraries
#LIBS = -L$(LAPACK_PATH) -llapack -L$(BLAS_PATH) -lblas
LIBS = -llapack -lblas

SRC = module_constants.f90 module_Numerical_recipes.f90 module_kernel.f90 module_LRE.f90  \
			module_set_up.f90 IO.f90 module_units.f90 module_constants.f90 module_reading_dep.f90 \
			module_input_output.f90 module_projection.f90 main_sphincs.f90


# Object files (replace .f90 with .o)
OBJ = $(SRC:.f90=.o)

# Executable name
EXEC = map_it

# Default target
all: $(EXEC)

# Link the objects to create the executable
$(EXEC): $(OBJ)
	$(FC) $(FLAGS) -o $(EXEC) $(OBJ) $(LIBS)

# Compile each source file into an object file
%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

# Clean the build directory
clean:
	rm -f $(OBJ) $(EXEC) *.mod

.PHONY: all clean


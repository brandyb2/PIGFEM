#
# Makefile for PCIGFEM
#

current_dir = $(shell pwd)


# PETSc includes and libraries
PETSC_VERSION		= 3.7.3
PETSC_DIR			= ${current_dir}/petsc-${PETSC_VERSION}
PETSC_ARCH_OPT		= arch_external_opt
PETSC_ARCH_DEBUG	= arch_external_debug
ifeq (${MAKECMDGOALS},debug)
	include ${PETSC_DIR}/${PETSC_ARCH_DEBUG}/lib/petsc/conf/petscvariables
else
	include ${PETSC_DIR}/${PETSC_ARCH_OPT}/lib/petsc/conf/petscvariables
endif
PETSC_LIB = ${PETSC_KSP_LIB}
PETSC_INCLUDE = ${PETSC_CC_INCLUDES}

# GSL includes and libraries
GSL_VERSION		= 2.1
GSL_INCLUDE		= -I${current_dir}/gsl-${GSL_VERSION}/install/include/gsl
GSL_LIB			= -Wl,-rpath -Wl,${current_dir}/gsl-${GSL_VERSION}/install/lib -lgsl -lgslcblas

# List of all includes and libraries
LOCAL_INCLUDES	:= -I${current_dir}/include
INCLUDES		:= ${LOCAL_INCLUDES} ${PETSC_CC_INCLUDES} ${GSL_INCLUDE}
LIBS			:= ${GSL_LIB} ${PETSC_LIB} -lm




# Basic compiler variable definitions
ifeq (${MAKECMDGOALS},debug)
	CXX			:= ${PETSC_DIR}/${PETSC_ARCH_DEBUG}/bin/mpicxx
	CXX_LINK	:= ${CXX}
else
	CXX			:= ${PETSC_DIR}/${PETSC_ARCH_OPT}/bin/mpicxx
	CXX_LINK	:= ${CXX}
endif
CXXFLAGS		:= -g -Wall -O3 -pg -std=c++0x
LDFLAGS			:= -O3 -pg
CXXFLAGS_DEBUG	:= -g -Wall -O0 -std=c++0x
LDFLAGS_DEBUG	:= -O0 -pg


# Define where all of the source code is coming from and where all of the object files are going
SOURCES			:= $(wildcard src/*.cpp)
OBJECTS			:= $(patsubst src/%.cpp,obj/%_o.o,$(SOURCES))
OBJECTS_DEBUG	:= $(patsubst src/%.cpp,obj/%_d.o,$(SOURCES))



# Define which execuutable I want to build
EXECUTABLE = Main
#EXECUTABLE = Main_FD_Shape
#EXECUTABLE = GoodierTest
#EXECUTABLE = CohesivePlanes

EXEC_SOURCE = Tests/${EXECUTABLE}.cpp
EXECUTABLE_DEBUG = ${EXECUTABLE}_debug
EXEC_OBJECT = obj/${EXECUTABLE}.o
EXEC_OBJECT_DEBUG = obj/${EXECUTABLE}_debug.o




# Default to building the optimized code
.PHONY: all
all: ${EXECUTABLE}

# Linker rule to link all object files to form the executable
${EXECUTABLE}: ${EXEC_OBJECT} ${OBJECTS}
	${CXX_LINK} ${LDFLAGS} -o $@ $^ ${LIBS}

# Special rule to form the EXEC object file because it comes from the tests directory
${EXEC_OBJECT}: ${EXEC_SOURCE}
	${CXX} ${CXXFLAGS} -c ${EXEC_SOURCE}  -o $@ ${INCLUDES}

# Pattern rule to form all object files for all .cpp files existing in the src directory
obj/%_o.o: src/%.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@ ${INCLUDES}




# Build the debug versions
debug: ${EXECUTABLE_DEBUG}

# Linker rule to link all object files to form the executable
${EXECUTABLE_DEBUG}: ${EXEC_OBJECT_DEBUG} ${OBJECTS_DEBUG}
	${CXX_LINK} ${LDFLAGS_DEBUG} -o $@ $^ ${LIBS}

# Special rule to form the EXEC object file because it comes from the tests directory
${EXEC_OBJECT_DEBUG}: ${EXEC_SOURCE}
	${CXX} ${CXXFLAGS_DEBUG} -c ${EXEC_SOURCE}  -o $@ ${INCLUDES}

# Pattern rule to form all object files for all .cpp files existing in the src directory
obj/%_d.o: src/%.cpp
	${CXX} ${CXXFLAGS_DEBUG} -c $< -o $@ ${INCLUDES}


check:
	${LIBS} \
	\
	\
	${INCLUDES}


# Clean the build directories
.PHONY: clean
clean:
	rm -f *.o *~       \
	rm -f src/*~       \
	rm -f include/*~   \
	rm -f Tests/*~     \
	rm -f obj/*.o

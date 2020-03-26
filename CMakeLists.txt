cmake_minimum_required(VERSION 3.0)

project(mpitest)



file(GLOB sources 
  "${PROJECT_SOURCE_DIR}/include/*.hpp",
  "${PROJECT_SOURCE_DIR}/src/*.cpp")




#################
## MKL ##########
#################

find_package(MKL REQUIRED)
find_package(Threads)

add_definitions(-DHAVE_MKL)
add_definitions(-DUSE_PARDISO)
set(HAVE_MKL 1)
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARY_DIR})

message("${MKL_INCLUDE_DIR} ----------------------------------")

include_directories(${METIS_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

set(BLAS_VENDOR Intel10_64lp) 
find_package(BLAS REQUIRED) 
find_package(LAPACK REQUIRED)

set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} m dl pthread mkl_core mkl_sequential)
set(MKL_LIB_DIR ${MKL_LIB_DIR} m dl mkl_sequential mkl_core pthread)

#################
## MPI ##########
#################
find_package(MPI REQUIRED)
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)

#################
## VTK ##########
#################
set(VTK_DIR  "$ENV{VTK_ROOT}/build" CACHE PATH "")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#################
## EIGEN ########
#################
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(EIGEN3_INCLUDE_DIR)

#################
## METIS ########
#################
set(METIS_LIBRARIES} $ENV{METIS_LIBRARIES} CACHE PATH "Where the MKL are stored")
set(METIS_INCLUDE_DIRS $ENV{METIS_INCLUDEDIRS} CACHE PATH "Where the MKL are stored")
include_directories(${METIS_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

# Use imported targets would make things much eazier. Thanks Levi for pointing it out.
#aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} main.cpp ${sources})
target_link_libraries(${PROJECT_NAME}  PRIVATE MPI::MPI_C  ${VTK_LIBRARIES} Eigen3::Eigen 
  ${METIS_LIBRARIES} )

set(CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "" FORCE)

set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

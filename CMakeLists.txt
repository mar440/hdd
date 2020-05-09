cmake_minimum_required(VERSION 3.0)

project(hdd)

add_compile_definitions(DBG=0 DMP_INPUTS=2)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/") 
file(GLOB sources 
  "${PROJECT_SOURCE_DIR}/include/*.hpp","${PROJECT_SOURCE_DIR}/*.hpp",
  "${PROJECT_SOURCE_DIR}/src/*.cpp")


if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -fPIC -Wunused-but-set-variable")
endif()


#################
## MKL ##########
#################
##
find_package(MKL REQUIRED ) 
find_package(Threads)

add_definitions(-DHAVE_MKL)
add_definitions(-DUSE_PARDISO)
set(HAVE_MKL 1)
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARY_DIR})

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
## EIGEN ########
#################
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(EIGEN3_INCLUDE_DIR)

#   #################
#   ## METIS ########
#   #################
#   set(METIS_LIBRARIES $ENV{METIS_LIBRARIES}/libmetis.a CACHE PATH "mydoc")
#   set(METIS_INCLUDEDIRS $ENV{METIS_INCLUDEDIRS} CACHE PATH "mydoc")
#   include_directories(${METIS_INCLUDEDIRS} ${CMAKE_CURRENT_SOURCE_DIR})

#################
## BOOST ########
#################

set(Boost_USE_STATIC_RUNTIME OFF) 

set(BOOST_PACKAGES "date_time" "filesystem" "system" "chrono" "thread" "regex" "unit_test_framework")

# Boost version: Windows32(1.55) and Linux64(1.52)
if(WIN32)
 if(VTK_MAJOR_VERSION GREATER 6)
  set(BOOST_VERSION_REQUIRED "1.64")
 else()
  set(BOOST_VERSION_REQUIRED "1.55")
 endif()

 add_definitions("/D _SCL_SECURE_NO_WARNINGS")
 #set(Boost_COMPILER -vc120)
 #add_definitions(-DBOOST_ALL_NO_LIB)
else(WIN32)
 set(BOOST_VERSION_REQUIRED "1.58")
 add_definitions(-DBOOST_LOG_DYN_LINK)
endif(WIN32)

add_definitions(-DBOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE)

if(WIN32)
  find_package(Boost ${BOOST_VERSION_REQUIRED} REQUIRED) # COMPONENTS ${BOOST_PACKAGES})
  link_directories(${Boost_LIBRARY_DIRS})
  link_directories(${BOOST_LIBRARYDIR})
  message(STATUS "BOOST LIBRARIES: ${Boost_LIBRARIES}")
  message(STATUS "BOOST LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
else(WIN32)
  find_package(Boost ${BOOST_VERSION_REQUIRED} COMPONENTS ${BOOST_PACKAGES} REQUIRED)
endif(WIN32)

add_definitions(-DBOOST_MINOR_VERSION=${Boost_MINOR_VERSION})

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
else()
  message("boost not!!! not found ----------------------------------")
endif()

set(CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON )

add_library(${PROJECT_NAME} SHARED main.cpp ${sources})
target_link_libraries(${PROJECT_NAME}  PRIVATE MPI::MPI_C   Eigen3::Eigen 
   ${BLAS_LIBRARIES}  ${Boost_LIBRARIES})
 #${METIS_LIBRARIES} 


install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
install(FILES hddApi.hpp hddPythonApi.hpp DESTINATION ${PROJECT_SOURCE_DIR})
#install(FILES hddPythonApi.hpp DESTINATION ${PROJECT_SOURCE_DIR})

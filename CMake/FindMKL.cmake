###     ################################################################################
###     #
###     # \file      cmake/FindMKL.cmake
###     # \author    J. Bakosi
###     # \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
###     # \brief     Find the Math Kernel Library from Intel
###     # \date      Thu 26 Jan 2017 02:05:50 PM MST
###     #
###     ################################################################################
###     
###     # Find the Math Kernel Library from Intel
###     #
###     #  MKL_FOUND - System has MKL
###     #  MKL_INCLUDE_DIRS - MKL include files directories
###     #  MKL_LIBRARIES - The MKL libraries
###     #  MKL_INTERFACE_LIBRARY - MKL interface library
###     #  MKL_SEQUENTIAL_LAYER_LIBRARY - MKL sequential layer library
###     #  MKL_CORE_LIBRARY - MKL core library
###     #
###     #  The environment variables MKLROOT and INTEL are used to find the library.
###     #  Everything else is ignored. If MKL is found "-DMKL_ILP64" is added to
###     #  CMAKE_C_FLAGS and CMAKE_CXX_FLAGS.
###     #
###     #  Example usage:
###     #
###     #  find_package(MKL)
###     #  if(MKL_FOUND)
###     #    target_link_libraries(TARGET ${MKL_LIBRARIES})
###     #  endif()
###     
###     # If already in cache, be silent
###     
###     
###     
###     if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES AND MKL_INTERFACE_LIBRARY AND
###         MKL_SEQUENTIAL_LAYER_LIBRARY AND MKL_CORE_LIBRARY)
###       set (MKL_FIND_QUIETLY TRUE)
###     endif()
###     
###     if(NOT BUILD_SHARED_LIBS)
###       set(INT_LIB "libmkl_intel_ilp64.a")
###       set(SEQ_LIB "libmkl_sequential.a")
###       set(THR_LIB "libmkl_intel_thread.a")
###       set(COR_LIB "libmkl_core.a")
###     else()
###       set(INT_LIB "mkl_intel_ilp64")
###       set(SEQ_LIB "mkl_sequential")
###       set(THR_LIB "mkl_intel_thread")
###       set(COR_LIB "mkl_core")
###     endif()
###     
###     find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS $ENV{MKLROOT}/include)
###     
###     find_library(MKL_INTERFACE_LIBRARY
###                  NAMES ${INT_LIB}
###                  PATHS $ENV{MKLROOT}/lib
###                        $ENV{MKLROOT}/lib/intel64
###                        $ENV{INTEL}/mkl/lib/intel64
###                  NO_DEFAULT_PATH)
###     
###     find_library(MKL_SEQUENTIAL_LAYER_LIBRARY
###                  NAMES ${SEQ_LIB}
###                  PATHS $ENV{MKLROOT}/lib
###                        $ENV{MKLROOT}/lib/intel64
###                        $ENV{INTEL}/mkl/lib/intel64
###                  NO_DEFAULT_PATH)
###     
###     find_library(MKL_CORE_LIBRARY
###                  NAMES ${COR_LIB}
###                  PATHS $ENV{MKLROOT}/lib
###                        $ENV{MKLROOT}/lib/intel64
###                        $ENV{INTEL}/mkl/lib/intel64
###                  NO_DEFAULT_PATH)
###     
###     set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
###     set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
###     
###     if (MKL_INCLUDE_DIR AND
###         MKL_INTERFACE_LIBRARY AND
###         MKL_SEQUENTIAL_LAYER_LIBRARY AND
###         MKL_CORE_LIBRARY)
###     
###         if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
###             NOT DEFINED ENV{CRAY_PRGENVGNU} AND
###             NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
###             NOT DEFINED ENV{CRAY_PRGENVINTEL})
###           set(ABI "-m64")
###         endif()
###     
###         set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64 ${ABI}")
###         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64 ${ABI}")
###     
###     else()
###     
###       set(MKL_INCLUDE_DIRS "")
###       set(MKL_LIBRARIES "")
###       set(MKL_INTERFACE_LIBRARY "")
###       set(MKL_SEQUENTIAL_LAYER_LIBRARY "")
###       set(MKL_CORE_LIBRARY "")
###     
###     endif()
###     
###     # Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
###     # all listed variables are TRUE.
###     INCLUDE(FindPackageHandleStandardArgs)
###     FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_INTERFACE_LIBRARY MKL_SEQUENTIAL_LAYER_LIBRARY MKL_CORE_LIBRARY)
###     
###     MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY MKL_SEQUENTIAL_LAYER_LIBRARY MKL_CORE_LIBRARY)
###     
###     


###     # - Find the MKL libraries
###     # Modified from Armadillo's ARMA_FindMKL.cmake
###     # This module defines
###     #  MKL_INCLUDE_DIR, the directory for the MKL headers
###     #  MKL_LIB_DIR, the directory for the MKL library files
###     #  MKL_COMPILER_LIB_DIR, the directory for the MKL compiler library files
###     #  MKL_LIBRARIES, the libraries needed to use Intel's implementation of BLAS & LAPACK.
###     #  MKL_FOUND, If false, do not try to use MKL; if true, the macro definition USE_MKL is added.
###     
###     # Set the include path
###     # TODO: what if MKL is not installed in /opt/intel/mkl?
###     # try to find at /opt/intel/mkl
###     # in windows, try to find MKL at C:/Program Files (x86)/Intel/Composer XE/mkl
###     
###     #if ( WIN32 )
###     #  if(NOT DEFINED ENV{MKLROOT})
###     #    set(MKLROOT "C:/Program Files (x86)/Intel/Composer XE" CACHE PATH "Where the MKL are stored")
###     #  endif(NOT DEFINED ENV{MKLROOT}) 
###     #else ( WIN32 )
###     #    set(MKLROOT "/home/mar440/intel/" CACHE PATH "Where the MKL are stored")
###     #endif ( WIN32 )
###     
###     

if(DEFINED ENV{MKLROOT})
    set(MKLROOT $ENV{MKLROOT} CACHE PATH "Where the MKL are stored")
endif()


if (EXISTS ${MKLROOT}/mkl)
    SET(MKL_FOUND TRUE)
    message("MKL is found at ${MKLROOT}/mkl")
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( USE_MKL_64BIT On )
        if ( ARMADILLO_FOUND )
            if ( ARMADILLO_BLAS_LONG_LONG )
                set( USE_MKL_64BIT_LIB On )
                ADD_DEFINITIONS(-DMKL_ILP64)
                message("MKL is linked against ILP64 interface ... ")
            endif ( ARMADILLO_BLAS_LONG_LONG )
        endif ( ARMADILLO_FOUND )
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( USE_MKL_64BIT Off )
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8)
else (EXISTS ${MKLROOT}/mkl)
    SET(MKL_FOUND FALSE)
    message("MKL is NOT found ... ")
endif (EXISTS ${MKLROOT}/mkl)

if (MKL_FOUND)
    set(MKL_INCLUDE_DIR "${MKLROOT}/mkl/include")
    ADD_DEFINITIONS(-DUSE_MKL)
    if ( USE_MKL_64BIT )
        set(MKL_LIB_DIR "${MKLROOT}/mkl/lib/intel64")
        set(MKL_COMPILER_LIB_DIR "${MKLROOT}/compiler/lib/intel64")
		set(COMPILER_REDIST "${MKLROOT}/redist/intel64/compiler")
        set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT}/lib/intel64")
        if ( USE_MKL_64BIT_LIB )
                if (WIN32)
                    set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_ilp64)
                else (WIN32)
                    set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_ilp64)
                endif (WIN32)
        else ( USE_MKL_64BIT_LIB )
                if (WIN32)
                    set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64)
                else (WIN32)
                    set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64)
                endif (WIN32)
        endif ( USE_MKL_64BIT_LIB )
    else ( USE_MKL_64BIT )
        set(MKL_LIB_DIR "${MKLROOT}/mkl/lib/ia32")
        set(MKL_COMPILER_LIB_DIR "${MKLROOT}/compiler/lib/ia32")
        set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT}/lib/ia32")
        if ( WIN32 )
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_c)
        else ( WIN32 )
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel)
        endif ( WIN32 )
    endif ( USE_MKL_64BIT )

    if (WIN32)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_thread)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_core)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} libiomp5md)
    else (WIN32)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_gnu_thread)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_core)
    endif (WIN32) 
endif (MKL_FOUND)

IF (MKL_FOUND)
    IF (NOT MKL_FIND_QUIETLY)
        MESSAGE(STATUS "Found MKL libraries: ${MKL_LIBRARIES}")
        MESSAGE(STATUS "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}")
        MESSAGE(STATUS "MKL_LIB_DIR: ${MKL_LIB_DIR}")
        MESSAGE(STATUS "MKL_COMPILER_LIB_DIR: ${MKL_COMPILER_LIB_DIR}")
    ENDIF (NOT MKL_FIND_QUIETLY)

    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
ELSE (MKL_FOUND)
    IF (MKL_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find MKL libraries")
    ENDIF (MKL_FIND_REQUIRED)
ENDIF (MKL_FOUND)

# MARK_AS_ADVANCED(MKL_LIBRARY)

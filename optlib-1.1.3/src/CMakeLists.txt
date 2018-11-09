include_directories(${optlib_SOURCE_DIR}/src)
include_directories(${optlib_BINARY_DIR}/src)
link_directories(${optlib_BINARY_DIR}/src)

set (optlib_VERSION_MAJOR 1)
set (optlib_VERSION_MINOR 1.3)

configure_file (
    "${optlib_SOURCE_DIR}/src/ol_config.h.in"
    "${optlib_BINARY_DIR}/src/ol_config.h"
)

option(USE_MPI "Compile mpi version of library" OFF)
set(CMAKE_VERBOSE_MAKEFILE OFF)
if(USE_MPI)
    set(LIBNAME opt_mpi)
    set(EXECNAME optlib_driver_mpi)
    find_package(MPI REQUIRED)
    mark_as_advanced(CLEAR MPI_INCLUDE_PATH
        MPI_COMPILE_FLAGS MPI_LINK_FLAGS
        MPI_LIBRARY MPI_EXTRA_LIBRARY)
    SET(CMAKE_C_COMPILER mpicc)
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS}
        ${MPI_COMPILE_FLAGS})
    add_definitions(-DUSE_OPTMPI)
    include_directories(MPI_INCLUDE_PATH)
    add_executable(${EXECNAME} ol_driver_mpi.c func.c)
    target_link_libraries(${EXECNAME} ${MPI_LIBRARIES})
else(USE_MPI)
    set(EXECNAME optlib_driver)
    set(LIBNAME opt)
    add_executable(${EXECNAME} ol_driver_serial.c func.c)
endif(USE_MPI)

add_library(${LIBNAME} ol_fw_mpi.c ol_genetic_alg.c ol_interp.c ol_mem.c
    ol_powell.c ol_rand_tools.c ol_simanneal.c ol_transform.c optlib.c)
target_link_libraries(${EXECNAME} ${LIBNAME})
if(UNIX)
  target_link_libraries(${EXECNAME} m)
endif(UNIX)

install(FILES ol_fw_mpi.h ol_genetic_alg.h ol_interp.h ol_mem.h ol_powell.h
    ol_rand_tools.h ol_simanneal.h ol_transform.h optlib.h optlib_mpi.h
${optlib_BINARY_DIR}/src/ol_config.h
    DESTINATION include)
install(TARGETS ${LIBNAME} ${EXECNAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

cmake_minimum_required (VERSION 2.6)
project(optlib)
add_subdirectory(src)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "/usr/local/optlib" CACHE PATH "FOO install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

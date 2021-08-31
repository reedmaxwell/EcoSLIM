cmake_minimum_required(VERSION 3.14)

project(ecoslim LANGUAGES Fortran)

find_package(OpenMP)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(DEBUG_FLAGS "-ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(DEBUG_FLAGS  "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(DEBUG_FLAGS "-C")
endif()

#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FLAGS}")

set(SOURCE ran1.mod.f90 EcoSLIM.f90 pfb_read.f90 vtk_write.f90	vtk_write_points.f90)

add_executable(EcoSLIM.exe ${SOURCE})

if(OpenMP_Fortran_FOUND)
  target_link_libraries(EcoSLIM.exe PUBLIC OpenMP::OpenMP_Fortran)
endif()


set(source
  src/module_enkf.f90
)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmkl_rt")

add_library(enkf ${source})
target_link_directories(enkf PUBLIC ${LAPACK_LIBDIR})


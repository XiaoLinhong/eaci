set(sources
  src/module_hash.F90
  src/module_constant.f90
  src/module_structure.f90

  src/module_csvio.f90
  src/module_neighbour.f90
  src/module_routine.f90
)

add_library(task STATIC ${sources})
target_link_libraries(task datetime tool flogger datetime)

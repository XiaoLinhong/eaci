
add_subdirectory(lib/container)
add_subdirectory(lib/string)

set(source_files
  src/face_mod.F90
  src/log_mod.F90
  src/flogger.F90
)

add_library(flogger ${source_files})
target_link_libraries(flogger fortran_container fortran_string)


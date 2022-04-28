
set(sources
  src/string_mod.F90
  src/string_actions_mod.F90
  src/string_numerics_mod.F90
  src/string.F90
)

add_library(fortran_string ${sources})


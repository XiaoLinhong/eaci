
set(source_files
  src/linked_list_mod.F90
  src/hash_table_mod.F90
)

add_library(fortran_container STATIC ${source_files})

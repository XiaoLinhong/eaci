set(sources
  src/module_datetime.f90
  src/module_timedelta.f90
)

add_library(datetime STATIC ${sources})

cmake_minimum_required(VERSION 3.11) 

project(eaci LANGUAGES Fortran)

MESSAGE( STATUS "test?: ${TEST_FLAGS}")

# 编译选项
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  # 累积编译选项
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -Ofast -qopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -Ofast")

  if (CMAKE_BUILD_TYPE MATCHES Debug) #
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  endif ()

elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

  if (CMAKE_BUILD_TYPE MATCHES Debug) # 测试
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()

endif ()

# 打印编译器
MESSAGE( STATUS "compiler is: ${CMAKE_Fortran_COMPILER_ID}")
MESSAGE( STATUS "compiler flag: ${CMAKE_Fortran_FLAGS}")
MESSAGE( STATUS "compiler type: ${CMAKE_BUILD_TYPE}")

# module位置
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
MESSAGE( STATUS "module directory: ${CMAKE_Fortran_MODULE_DIRECTORY}")

# 获取父目录, 为了嵌套在其他的仓中进行编译
get_directory_property(parent_dir PARENT_DIRECTORY)
MESSAGE( STATUS "parent directory: ${parent_dir}")

# 编译库
add_subdirectory(lib/enkf)
add_subdirectory(lib/tool)
add_subdirectory(lib/flogger)
add_subdirectory(lib/datetime)

add_subdirectory(lib/task)

set(sources
    src/module_cfg.f90
    src/main.f90
)

# 添加数学库
set(LAPACK_LIBDIR "/public/home/bedrock/envs/v1.0/intel/2019.4/mkl/lib/intel64_lin/lib")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmkl_rt")

add_executable(eaci.exe ${sources})
target_link_directories(eaci.exe PUBLIC ${LAPACK_LIBDIR})
target_link_libraries(eaci.exe task tool flogger datetime enkf)

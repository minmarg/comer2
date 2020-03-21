# refer to the root source directory of the project as ${<PRJNAME>_SOURCE_DIR} and
# to the root binary directory of the project as ${<PRJNAME>_BINARY_DIR}.
cmake_minimum_required (VERSION 3.8)
project (COMER2 LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

## set build configuration
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

if(WIN32)
  add_compile_definitions(_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES=1)
endif()

## set definitions for configuration DEBUG, required when this 
## configuration is in use
#set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_DEBUG "__DEBUG__")
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS $<$<CONFIG:Debug>:__DEBUG__>)


if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set_property(DIRECTORY APPEND PROPERTY COMPILE_OPTIONS $<$<CONFIG:Debug>:-O2>)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -Wall -Xcompiler -Wextra -Xptxas -v")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  ##set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  ##add_compile_options(-Wall -Wextra)
endif()



FOREACH(LANG C CXX CUDA CUDACXX)
  ##set(CMAKE_${LANG}_FLAGS_DEBUG 
  ##  "${CMAKE_${LANG}_FLAGS_DEBUG} -D__DEBUG__")
  ## source file extension is replaced when producing object files
  set(CMAKE_${LANG}_OUTPUT_EXTENSION_REPLACE 1)
ENDFOREACH()

## set local directories (bin,state)
set(LOCALBINDIR bin)
set(LOCALSTATEDIR var)

# Recurse into the subdirectories. 
# The same process will walk through the project's entire directory structure.
## libs...
add_subdirectory (extsp)
add_subdirectory (liblib)
add_subdirectory (libseg)
add_subdirectory (libHDP)
add_subdirectory (libpro)
add_subdirectory (libmycu)
## execs...
add_subdirectory (makepro_and_run_comer)
add_subdirectory (comer2msa)
add_subdirectory (blast2fa)
add_subdirectory (makepro)
add_subdirectory (makecov)
add_subdirectory (makedb)
add_subdirectory (db2bin)
add_subdirectory (comer)
## configuration...
add_subdirectory (env)


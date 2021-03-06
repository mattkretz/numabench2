cmake_minimum_required(VERSION 2.8.3)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

project(numabench2)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebug RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(Vc 0.7.0 REQUIRED)

# We want C++11 code
AddCompilerFlag("-std=c++11" _yes)
if(NOT _yes)
   AddCompilerFlag("-std=c++0x")
endif()

# need pthreads
AddCompilerFlag(-pthread)

# disable automatic addition of software prefetch instructions
AddCompilerFlag("-fno-prefetch-loop-arrays")

find_file(HAVE_SYS_MMAN sys/mman.h)
if(NOT HAVE_SYS_MMAN)
   message(ERROR "could not find required sys/mman.h header to enable memlocking")
   add_definitions(-DVC_BENCHMARK_NO_MLOCK)
endif()
mark_as_advanced(HAVE_SYS_MMAN)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   set(REALTIME_BENCHMARKS FALSE CACHE BOOL "Use SCHED_FIFO priority to give benchmarks highest priority, achieving less jitter in the results.")
   if(REALTIME_BENCHMARKS)
      add_executable(rtwrapper rtwrapper.cpp)
      add_definitions(-DSCHED_FIFO_BENCHMARKS)
      add_custom_command(TARGET rtwrapper POST_BUILD
         COMMAND sudo chown root rtwrapper
         COMMAND sudo chmod +s rtwrapper
         WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
         COMMENT "Making the rtwrapper helper suid root"
         )
   endif(REALTIME_BENCHMARKS)
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

include_directories(${Vc_INCLUDE_DIR})
add_definitions(${Vc_DEFINITIONS})

add_executable(cpuid cpuid.cpp)
target_link_libraries(cpuid ${Vc_LIBRARIES})

add_library(cpuset STATIC cpuset.c)

add_executable(add1test add1test.cpp)
target_link_libraries(add1test ${Vc_LIBRARIES} rt cpuset)

add_executable(numabench numabench.cpp)
target_link_libraries(numabench ${Vc_LIBRARIES} rt cpuset)

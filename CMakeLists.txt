cmake_minimum_required(VERSION 2.8.3)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BUILD_DIR)
   message(ERROR "You don't want to configure in the source directory!")
endif()

project(numabench2)

find_package(Vc REQUIRED)

# We want C++11 code
AddCompilerFlag("-std=c++11" _yes)
if(NOT _yes)
   AddCompilerFlag("-std=c++0x")
endif()

find_file(HAVE_SYS_MMAN sys/mman.h)
if(NOT HAVE_SYS_MMAN)
   message(ERROR "could not find required sys/mman.h header to enable memlocking")
   add_definitions(-DVC_BENCHMARK_NO_MLOCK)
endif()
mark_as_advanced(HAVE_SYS_MMAN)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   set(REALTIME_BENCHMARKS FALSE CACHE BOOL "Use SCHED_FIFO priority to give benchmarks highest priority, achieving less jitter in the results.")
   if(REALTIME_BENCHMARKS)
      add_executable(benchmark benchmark.cpp)
      add_definitions(-DSCHED_FIFO_BENCHMARKS)
      add_custom_command(TARGET benchmark POST_BUILD
         COMMAND sudo chown root benchmark
         COMMAND sudo chmod +s benchmark
         WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
         COMMENT "Making the benchmark helper suid root"
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

cmake_minimum_required(VERSION 3.13)

project(colvars_sde CXX)

get_filename_component(COLVARS_SDE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

option(COLVARS_SDE_DEBUG "Build debugging binaries" OFF)

file(GLOB COLVARS_SDE_SOURCES ${COLVARS_SDE_SOURCE_DIR}/src/[^.]*.cpp)

add_executable(colvars_sde ${COLVARS_SDE_SOURCES})

option(COLVARS_TORCH "Link against the Torch library" OFF)

if(COLVARS_TORCH)
  if (CMAKE_CXX_STANDARD LESS 17)
    message(STATUS "Increase the CMAKE_CXX_STANDARD to 17")
    set_property(TARGET colvars_sde PROPERTY CXX_STANDARD 17)
  else()
    message(STATUS "Keep the CMAKE_CXX_STANDARD as ${CMAKE_CXX_STANDARD}")
  endif()
  if(DEFINED LIBTORCH_PREFIX)
    find_package(Torch REQUIRED PATHS file(JOIN ${LIBTORCH_PREFIX} "share/cmake/Torch"))
  else()
    find_package(Torch REQUIRED)
  endif()
  if(Torch_FOUND)
    target_compile_definitions(colvars_sde PRIVATE -DCOLVARS_TORCH)
    target_compile_options(colvars_sde PRIVATE ${TORCH_CXX_FLAGS})
    target_include_directories(colvars_sde PRIVATE ${TORCH_INCLUDE_DIRS})
    target_link_libraries(colvars_sde PUBLIC ${TORCH_LIBRARIES})
  endif()
endif()


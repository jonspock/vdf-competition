
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

project(VDF)

if(${CMAKE_SYSTEM_NAME} MATCHES "XXXDarwin")
  set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/7.0.0/bin/clang")
  set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/7.0.0/bin/clang++")
  set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/7.0.0/lib")
  set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/7.0.0/include")
endif()


# Add path for custom modules
set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH}
	${CMAKE_CURRENT_SOURCE_DIR}/modules
)

enable_testing()

find_package(gmp)
if (GMP_FOUND)
message(STATUS "Found libgmp")
endif()

find_package(sodium REQUIRED)
if (SODIUM_FOUND)
  message(STATUS "Found libsodium")
  message(STATUS "Sodium include dir = ${sodium_INCLUDE_DIR}")
endif()


include(CheckSymbolExists)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++17 HAVE_FLAG_STD_CXX17)
if (HAVE_FLAG_STD_CXX17)
  set(CMAKE_CXX_STANDARD 17)
  message(STATUS "Compiler supports C++17")
else()
  message(ERROR "Compiler must at least support C++17")
endif()


if (MP_CONFIG_OK)
  add_definitions(-fopenmp=libomp)
  set(OpenMP_C_FLAGS "-fopenmp=libomp")
  set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set(OpenMP_CXX_FLAGS "-fopenmp=libomp")
  set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  
  find_package(OpenMP)
  if (FOUND_OPENMP)
    message(STATUS "Found openmp")
  else()
    message(STATUS " openmp not found!")
  endif()
endif()

include_directories(
  # for omp.h
  #  /usr/local/Cellar/llvm/7.0.0/lib/clang/7.0.0/include

  ${GMP_INCLUDE_DIR}
  ${sodium_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src/relic/include
  ${CMAKE_CURRENT_SOURCE_DIR}
  # For relic_conf.h file
  ${CMAKE_CURRENT_BINARY_DIR}
  )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/..)

#add_definitions(-mtune=haswell -march=haswell)

add_subdirectory(src)


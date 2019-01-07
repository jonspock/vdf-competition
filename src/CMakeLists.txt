CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RELEASE")
ENDIF()

#------------------- ------------------- ------------------- ------------------- -------------------
# For Relic
#------------------- ------------------- ------------------- ------------------- -------------------

##SET(CMAKE_VERBOSE_MAKEFILE OFF)
set(WORD 64 CACHE INTEGER "")
set(TIMER "CYCLE" CACHE STRING "")
set(CHECK "off" CACHE STRING "")
set(VERBS "off" CACHE STRING "")
set(SHLIB "OFF" CACHE STRING "")

set(FP_PRIME 381 CACHE INTEGER "")

IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(DSEED "UDEV" CACHE STRING "")
  set(STBIN "OFF" CACHE STRING "")
  set(FP_QNRES "off" CACHE STRING "")
ELSE()
  set(DSEED "DEV" CACHE STRING "")
  set(STBIN "OFF" CACHE STRING "")
  set(FP_QNRES "on" CACHE STRING "")
ENDIF()

# FOR GCD STEIN/LEHME seem good
# MXP MONTY/SLIDE same > BASIC
# MOD doesn't seem to matter -> use MONTY
set(BN_METHD "COMBA;COMBA;MONTY;SLIDE;LEHME;BASIC")
#set(BN_KARAT 4 CACHE INTERGER "") #-> Makes it worse

SET(BN_PRECI 2048)

set(FP_METHD "INTEG;INTEG;INTEG;MONTY;LOWER;SLIDE" CACHE STRING "")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(COMP "-g -O0 -funroll-loops -fomit-frame-pointer" CACHE STRING "")
else()
  set(COMP "-O3 -funroll-loops -fomit-frame-pointer" CACHE STRING "")
endif()


# Disable relic tests and benchmarks
set(TESTS 0 CACHE INTEGER "")
set(BENCH 0 CACHE INTEGER "")


file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
source_group("Headers" FILES ${HEADERS})

#add_definitions(-DARITH=GMP)

set(SOURCES
  classgroupSquarer.cpp
  proof_wesolowski.cpp
  )

add_library(vdf ${HEADERS} ${SOURCES})

add_subdirectory(relic)

set(tests
  cg
  bi
  bm
  pw
  mul
  compare
  sample
  sqr
  gcdext
  reduce
  mul_thread
  ###  mul_openmp
  mul_smalla
  vdf4
  dutta
  profile
  #  quick
  )

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  LINK_DIRECTORIES(/usr/local/Cellar/llvm/7.0.0/lib)
endif()

FOREACH(test ${tests})
	SET(SRC test_${test}.cpp)
	ADD_EXECUTABLE(test_${test} ${SRC})
  TARGET_LINK_LIBRARIES(test_${test} vdf sodium relic_s ${GMP_LIBRARIES} pthread)
  #set_property(TARGET test_${test} PROPERTY POSITION_INDEPENDENT_CODE FALSE)
ENDFOREACH(test ${tests})

#add_test(NAME benchmark COMMAND benchmark)

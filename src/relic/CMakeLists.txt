cmake_minimum_required(VERSION 3.1)
if(NOT ${CMAKE_VERSION} VERSION_LESS "3.1")
	cmake_policy(SET CMP0054 NEW)
endif()

project(RELIC C CXX)

set(PROJECT_VERSION_MAJOR "0")
set(PROJECT_VERSION_MINOR "4")
set(PROJECT_VERSION_PATCH "1")
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(VERSION ${PROJECT_VERSION})

set(INCLUDE ${CMAKE_CURRENT_BINARY_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include/low)
include_directories(${INCLUDE})

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

set(WFLAGS "-Wall")

message("-- Configuring ${PROJECT_NAME} ${PROJECT_VERSION}...")

  
option(DEBUG "Build with debugging support" off)
option(PROFL "Build with profiling support" off)
option(CHECK "Build with error-checking support" on)
option(VERBS "Build with detailed error messages" on)
option(TRACE "Build with tracing support" off)
option(OVERH "Build with overhead estimation" off)
option(DOCUM "Build documentation" on)
option(STRIP "Build only the selected algorithms" off)
option(QUIET "Build with printing disabled" off)
option(COLOR "Build with colored output" on)
option(BIGED "Build with big-endian support" off)
option(SHLIB "Build shared library" on)
option(STLIB "Build static library" on)
option(STBIN "Build static binaries" on)
option(AMALG "Build amalgamation" off)

if (SHOW_MESSAGES) 
  message(STATUS "Available switches (default = CHECK, VERBS, DOCUM):")
  
  message("   DEBUG=[off|on] Build with debugging support.")
  message("   PROFL=[off|on] Build with profiling support.")
  message("   CHECK=[off|on] Build with error-checking support.")
  message("   VERBS=[off|on] Build with detailed error messages.")
  message("   TRACE=[off|on] Build with tracing support.")
  message("   OVERH=[off|on] Build with overhead estimation.")
  message("   DOCUM=[off|on] Build documentation.")
  message("   STRIP=[off|on] Build only selected algorithms.")
  message("   QUIET=[off|on] Build with printing disabled.")
  message("   COLOR=[off|on] Build with colored output.")
  message("   BIGED=[off|on] Build with big-endian support.")
  message("   SHLIB=[off|on] Build shared library.")
  message("   STLIB=[off|on] Build static library.")
  message("   STBIN=[off|on] Build static binaries.")
  message("   AMALG=[off|on] Build amalgamation for better performance.")

  message(STATUS "Number of available processor cores (default = 1):")
  
  message("   CORES=n        If n > 1, please enable multithreading support.")
  
  message(STATUS "Available modules (default = ALL)")
  
  message("   WITH=BN       Multiple precision arithmetic.")
  message("   WITH=MD       Message digests (hash functions).")
  message("   Note: the programmer is responsible for not using unselected modules.")
  
  message(STATUS "Available arithmetic backends (default = easy):")

  message("   ARITH=easy     Easy-to-understand implementation.")
  message("   ARITH=gmp      GNU Multiple Precision library.\n")
  
  message(STATUS "Supported operating systems (default = LINUX):")
  
  message("   OPSYS=         Undefined/No specific operating system.")
  message("   OPSYS=LINUX    GNU/Linux operating system.")
  message("   OPSYS=FREEBSD  FreeBSD operating system.")
  message("   OPSYS=MACOSX   Mac OS X operating system.")
  message("   OPSYS=WINDOWS  Windows operating system.")
  message("   OPSYS=DROID    Android operating system.")
  message("   OPSYS=DUINO    Arduino platform.\n")
  
  message(STATUS "Supported multithreading APIs (default = UNDEF):")
  
  message("   MULTI=         No multithreading support.")
  message("   MULTI=OPENMP   Open Multi-Processing.")
  message("   MULTI=PTHREAD  POSIX threads.\n")

  message(STATUS "Prefix to identify this build of the library (default = \"\"):")
  
  message("   LABEL=relic\n")
endif()


include(cmake/arch.cmake)
include(cmake/bn.cmake)

# Number of cores.
set(CORES 1 CACHE INTEGER "Number of available processor cores.")

# Default modules.
set(WITH "ALL" CACHE STRING "Selected modules")
LIST(FIND WITH "ALL" TEMP)
if(TEMP GREATER -1)
	set(WITH_BN 1)
	set(WITH_MD 1)
endif(TEMP GREATER -1)

# Check if multiple precision integer arithmetic is required.
LIST(FIND WITH "BN" TEMP)
if(TEMP GREATER -1)
	set(WITH_BN 1)
endif(TEMP GREATER -1)

# Check if support for hash functions is required.
LIST(FIND WITH "MD" TEMP)
if(TEMP GREATER -1)
	set(WITH_MD 1)
endif(TEMP GREATER -1)

# Choose the arithmetic backend.
set(ARITH "gmp" CACHE STRING "Arithmetic backend")

# Compiler flags.
if(COMP)
  message(STATUS "COMP value was set at ${COMP}")
  # Not sure why CACHE used here, disable
else()
  set(COMP "-O3 -Wno-unused-function -funroll-loops -fomit-frame-pointer")
  message(STATUS "COMP value was not set, setting to ${COMP}")
endif()

# Linker flags.
string(TOLOWER ${ARITH} LFLAGS)
set(LFLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/src/low/${LFLAGS}/")
set(LINK "$ENV{LINK}" CACHE STRING "User-chosen linker flags.")
set(LFLAGS "${LFLAGS} ${LINK}")

if(STBIN)
	set(LFLAGS "${LFLAGS} -static")
endif(STBIN)

# Discover the operating system.
if(CMAKE_SYSTEM_NAME STREQUAL Linux)
	set(OPSYS "LINUX" CACHE STRING "Operating system")
else(CMAKE_SYSTEM_NAME STREQUAL Linux)
	if(CMAKE_SYSTEM_NAME STREQUAL FreeBSD)
		set(OPSYS "FREEBSD" CACHE STRING "Operating system")
	else(CMAKE_SYSTEM_NAME STREQUAL FreeBSD)
		if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
			set(OPSYS "MACOSX" CACHE STRING "Operating system")
		else(CMAKE_SYSTEM_NAME STREQUAL Darwin)
			set(OPSYS "RELIC_NONE" CACHE STRING "Operating system")
		endif(CMAKE_SYSTEM_NAME STREQUAL Darwin)
	endif(CMAKE_SYSTEM_NAME STREQUAL FreeBSD)
endif(CMAKE_SYSTEM_NAME STREQUAL Linux)
if(OPSYS STREQUAL LINUX)
	add_definitions(-D_GNU_SOURCE)
endif(OPSYS STREQUAL LINUX)
message(STATUS "Configured operating system: ${OPSYS}")

message(STATUS "COMP =  ${COMP}")

if(DEBUG)
	# If the user did not specify compile flags, we turn off all optimizations.
	set(CFLAGS "-O0 -fno-omit-frame-pointer")
	set(DFLAGS "-ggdb")
else(DEBUG)
	# If the user did not specify compile flags, we use sane defaults.
	set(CFLAGS "${COMP}")
	set(DFLAGS "")
endif(DEBUG)


if(MULTI STREQUAL OPENMP)
	find_package(OpenMP REQUIRED)
	set(CFLAGS "${CFLAGS} -fopenmp")
	set(MULTI "OPENMP" CACHE STRING "Multithreading interface")
else(MULTI STREQUAL OPENMP)
	if(MULTI STREQUAL PTHREAD)
		find_package(Threads REQUIRED)
		set(CFLAGS "${CFLAGS} -pthread")
		set(MULTI "PTHREAD" CACHE STRING "Multithreading interface")
	else(MULTI STREQUAL PTHREAD)
		set(MULTI "RELIC_NONE" CACHE STRING "Operating system")
	endif(MULTI STREQUAL PTHREAD)
endif(MULTI STREQUAL OPENMP)

if(ARITH STREQUAL "gmp")
	if(GMP_FOUND)
		include_directories(${GMP_INCLUDE_DIR})
		set(ARITH_LIBS ${GMP_LIBRARIES})
	endif(GMP_FOUND)
endif(ARITH STREQUAL "gmp")

set(CMAKE_C_FLAGS "-pipe -std=c99 ${AFLAGS} ${WFLAGS} ${DFLAGS} ${PFLAGS} ${CFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS ${LFLAGS})

message(STATUS "Compiler flags: ${CMAKE_C_FLAGS}")
message(STATUS "Linker flags: ${LFLAGS}")

string(TOUPPER ${ARITH} ARITH)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/relic_conf.h.in ${CMAKE_BINARY_DIR}/relic_conf.h @ONLY)
message(STATUS "Configured ${CMAKE_CURRENT_SOURCE_DIR}/include/relic_conf.h.in")
string(TOLOWER ${ARITH} ARITH)

if (LABEL)
	set(RELIC "relic_${LABEL}")
	set(RELIC_S "relic_s_${LABEL}")
else(LABEL)
	set(RELIC "relic")
	set(RELIC_S "relic_s")
endif(LABEL)

# Choose the arithmetic backend.
set(LABEL "" CACHE STRING "Build label")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB includes "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
install(FILES ${includes} DESTINATION include/${RELIC})
file(GLOB includes "${CMAKE_CURRENT_SOURCE_DIR}/include/low/*.h")
install(FILES ${includes} DESTINATION include/${RELIC}/low)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/" DESTINATION include/${RELIC})

install(FILES cmake/relic-config.cmake DESTINATION cmake/)

if(DOCUM)
	include(cmake/doxygen.cmake)
endif(DOCUM)


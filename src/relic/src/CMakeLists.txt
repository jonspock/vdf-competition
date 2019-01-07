string(TOLOWER ${ARITH} ARITH_PATH)
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/low/${ARITH_PATH}/CMakeLists.txt")
	include("${CMAKE_CURRENT_SOURCE_DIR}/low/${ARITH_PATH}/CMakeLists.txt")
endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/low/${ARITH_PATH}/CMakeLists.txt")

set(CORE_SRCS relic_err.c relic_core.c relic_util.c)

if (ARCH)
	string(TOLOWER ${ARCH} ARCH_PATH)
	if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/arch/relic_arch_${ARCH_PATH}.c")
		list(APPEND CORE_SRCS "arch/relic_arch_${ARCH_PATH}.c")
	endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/arch/relic_arch_${ARCH_PATH}.c")
else(ARCH)
	list(APPEND CORE_SRCS "arch/relic_arch_none.c")
endif(ARCH)

set(BN_SRCS
  bn/relic_bn.c
  )

set(MD_SRCS
  md/relic_md_sha256.c
  md/sha224-256.c)

set(ARITH_SRCS
  ./low/gmp/relic_bn_add_low.c
  ./low/gmp/relic_bn_cmp_low.c
  ./low/gmp/relic_bn_div_low.c
  ./low/gmp/relic_bn_mod_low.c
  ./low/gmp/relic_bn_mul_low.c
  ./low/gmp/relic_bn_shift_low.c
  ./low/gmp/relic_bn_sqr_low.c
  )

set(RELIC_SRCS ${CORE_SRCS} ${BN_SRCS} ${MD_SRCS} ${ARITH_SRCS})

# Not used
list(APPEND RELIC_SRCS ${ARITH_ASMS})
set_source_files_properties(${ARITH_ASMS} PROPERTIES LANGUAGE C COMPILE_FLAGS "-DASM -x assembler-with-cpp")

add_custom_target(arith_objs DEPENDS ${ARITH_OBJS})

macro(LINK_LIBS LIBRARY)
	if(OPSYS STREQUAL LINUX)
		target_link_libraries(${LIBRARY} rt)
	endif(OPSYS STREQUAL LINUX)
	if(MULTI STREQUAL OPENMP)
		target_link_libraries(${LIBRARY} gomp m)
	endif(MULTI STREQUAL OPENMP)
	if(MULTI STREQUAL PTHREAD)
		target_link_libraries(${LIBRARY} pthread m)
	endif(MULTI STREQUAL PTHREAD)
endmacro(LINK_LIBS)

function(AMALGAMATE output inputFiles)
	file(WRITE "${output}" "// This file is an automatic amalgamation of all relic source files.\n" )
	FOREACH(file ${inputFiles})
		file(APPEND "${output}" "#include \"${file}\"\n")
	ENDFOREACH()
endfunction(AMALGAMATE)

if(AMALG)
	AMALGAMATE("relic-all.c" "${RELIC_SRCS}${ARITH_OBJS}")
endif(AMALG)

if (STLIB)
	if(AMALG)
		add_library(${RELIC_S} STATIC "relic-all.c")
	else(AMALG)
		add_library(${RELIC_S} STATIC ${RELIC_SRCS} ${ARITH_OBJS})
	endif(AMALG)
	add_dependencies(${RELIC_S} arith_objs)
	target_link_libraries(${RELIC_S} ${ARITH_LIBS})
	LINK_LIBS(${RELIC_S})
	install(TARGETS ${RELIC_S} ARCHIVE DESTINATION lib)
endif(STLIB)

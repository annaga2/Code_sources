cmake_minimum_required (VERSION 2.6)
project (UTILZREG)


# Mandatory???
include_directories ("${PROJECT_BINARY_DIR}")
# end mandatory???

include_directories ("${PROJECT_SOURCE_DIR}/lib/nifti/niftilib")
include_directories ("${PROJECT_SOURCE_DIR}/lib/nifti/znzlib")
include_directories ("${PROJECT_SOURCE_DIR}/lib/include")
add_subdirectory (lib/nifti)
add_subdirectory (lib/src)

set (EXTRA_LIBS ${EXTRA_LIBS} Demons)
set (EXTRA_LIBS ${EXTRA_LIBS} Demons_plain)
set (EXTRA_LIBS ${EXTRA_LIBS} LDDMM)
set (EXTRA_LIBS ${EXTRA_LIBS} LDDMM_plain)
set (EXTRA_LIBS ${EXTRA_LIBS} GeoShoot)
set (EXTRA_LIBS ${EXTRA_LIBS} LIDM)
set (EXTRA_LIBS ${EXTRA_LIBS} SciCalcPack) #  SciCalcPack must be added after the other libs in EXTRA_LIBS.


IF(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")  #FOR OPENMP
ENDIF(NOT APPLE)


# add the executable
add_executable (uTIlzReg_Demons           applications/uTIlzReg_Demons.cc)
add_executable (uTIlzReg_Demons_plain     applications/uTIlzReg_Demons_plain.cc)
add_executable (uTIlzReg_LDDMM            applications/uTIlzReg_LDDMM.cc)
add_executable (uTIlzReg_LDDMM_plain      applications/uTIlzReg_LDDMM_plain.cc)
add_executable (uTIlzReg_GeoShoot         applications/uTIlzReg_GeoShoot.cc)
add_executable (uTIlzReg_LIDM             applications/uTIlzReg_LIDM.cc)
add_executable (uTIlzReg_Tools            applications/uTIlzReg_Tools.cc)

target_link_libraries (uTIlzReg_Demons  ${EXTRA_LIBS})
target_link_libraries (uTIlzReg_Demons_plain  ${EXTRA_LIBS})
target_link_libraries (uTIlzReg_LDDMM    ${EXTRA_LIBS})
target_link_libraries (uTIlzReg_LDDMM_plain  ${EXTRA_LIBS})
target_link_libraries (uTIlzReg_GeoShoot  ${EXTRA_LIBS})
target_link_libraries (uTIlzReg_LIDM  ${EXTRA_LIBS})
target_link_libraries (uTIlzReg_Tools  ${EXTRA_LIBS})

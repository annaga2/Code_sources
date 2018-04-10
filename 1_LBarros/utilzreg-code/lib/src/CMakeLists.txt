IF(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")  #FOR OPENMP
ENDIF(NOT APPLE)

add_library (Demons           Demons.cc) 
add_library (Demons_plain     Demons_plain.cc) 
add_library (LDDMM            LDDMM.cc)
add_library (LDDMM_plain      LDDMM_plain.cc) 
add_library (GeoShoot         GeoShoot.cc)
add_library (LIDM             LIDM.cc)
add_library (SciCalcPack      SciCalcPack.cc)

PROJECT(IJKMCUBE)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED IJK_DIR)
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED IJK_DIR)

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMCUBE_DIR "src/ijkmcube")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")
SET(IJK_ISOTABLE_DIR "${IJK_DIR}/isotable" CACHE PATH "Isotable directory")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(expat NrrdIO z)
ADD_DEFINITIONS(-DIJK_ISOTABLE_DIR=\"${IJK_ISOTABLE_DIR}\")

# Possibly turn on usage of zlib compression (requires linking with libz)
# (i.e., programs compiled with ITKNrrdIO must also be compiled with zlib)
ADD_DEFINITIONS(-DTEEM_ZLIB=1)

# Turn on TEEM_BUILD so that the proper dll export def's are 
# used on windows builds.
ADD_DEFINITIONS(-DTEEM_BUILD=1)

SET(IJKMCUBE_LIB_SRC ijkmcubeIO.cpp ijkmcube.cpp 
                     ijkmcube_datastruct.cpp ijkmcube_sub.cpp 
                     ijkmcube_extract.cpp ijkmcube_util.cpp ijksnapmc.cpp 
                     ijkMCtable.cpp ijkMCtable_poly.cpp ijkMCtable_ambig.cpp
                     ijkoctree.cpp ijkMCtable_xitIO.cpp)

SET(CMAKE_CXX_FLAGS "-std=c++20")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

ADD_EXECUTABLE(ijkmcube ijkmcube_main.cpp ${IJKMCUBE_LIB_SRC} )
  
ADD_LIBRARY(ijkmcubeL STATIC EXCLUDE_FROM_ALL ${IJKMCUBE_LIB_SRC} )
  
SET_TARGET_PROPERTIES(ijkmcubeL PROPERTIES OUTPUT_NAME ijkmcube)
ADD_CUSTOM_TARGET(lib DEPENDS ijkmcubeL)


SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkmcube DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkmcube DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkmcube DESTINATION "bin")
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()


ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKMCUBE_DIR}/ijkmcube.tar ${IJKMCUBE_DIR}/README ${IJKMCUBE_DIR}/INSTALL ${IJKMCUBE_DIR}/RELEASE_NOTES ${IJKMCUBE_DIR}/*.cpp ${IJKMCUBE_DIR}/*.h ${IJKMCUBE_DIR}/*.tpp ${IJKMCUBE_DIR}/CMakeLists.txt ${IJKMCUBE_DIR}/man/* ${IJKMCUBE_DIR}/ijkmcube_doxygen.config)

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkmcube_doxygen.config)

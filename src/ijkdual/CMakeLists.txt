PROJECT(IJKDUAL)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED IJK_DIR)
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED IJK_DIR)

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKDUAL_DIR "src/ijkdual")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include" "${IJK_DIR}/src/ijkdualtable")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(NrrdIO z)

# Possibly turn on usage of zlib compression (requires linking with libz)
# (i.e., programs compiled with ITKNrrdIO must also be compiled with zlib)
ADD_DEFINITIONS(-DTEEM_ZLIB=1)

# Turn on TEEM_BUILD so that the proper dll export def's are 
# used on windows builds.
ADD_DEFINITIONS(-DTEEM_BUILD=1)

ADD_COMPILE_OPTIONS(-Wall -Wno-unused-function)

# Include random library.
SET(CMAKE_CXX_FLAGS "-std=c++17 -D INCLUDE_RANDOM")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")


ADD_EXECUTABLE(ijkdual ijkdual_main.cpp ijkdualIO.cpp
                       ijkdual_datastruct.cpp)
  

SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkdual DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkdual DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkdual DESTINATION "bin")
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()


ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKDUAL_DIR}/ijkdual.tar ${IJKDUAL_DIR}/*.cpp ${IJKDUAL_DIR}/*.h ${IJKDUAL_DIR}/*.tpp ${IJKDUAL_DIR}/CMakeLists.txt)

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkdual_doxygen.config)

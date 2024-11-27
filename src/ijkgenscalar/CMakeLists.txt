PROJECT(IJKGENSCALAR)

#---------------------------------------------------------

SET(SRC_DIRNAME "src/ijkgenscalar")

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(NRRD_LIBDIR "${IJK_DIR}/lib")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(NrrdIO z)

# Possibly turn on usage of zlib compression (requires linking with libz)
# (i.e., programs compiled with ITKNrrdIO must also be compiled with zlib)
ADD_DEFINITIONS(-DTEEM_ZLIB=1)

# Turn on TEEM_BUILD so that the proper dll export def's are 
# used on windows builds.
ADD_DEFINITIONS(-DTEEM_BUILD=1)

SET(CMAKE_CXX_FLAGS "-std=c++20")

ADD_EXECUTABLE(ijkgenscalar ijkgenscalar.cpp ijkgenscalarIO.cpp)

SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkgenscalar DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkgenscalar DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkgenscalar DESTINATION "bin")
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()

# OBSOLETE
# SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
# INSTALL(TARGETS ijkgenscalar DESTINATION "bin")

SET(program_files ${SRC_DIRNAME}/README ${SRC_DIRNAME}/RELEASE_NOTES
                  ${SRC_DIRNAME}/*.cxx ${SRC_DIRNAME}/*.h 
                  ${SRC_DIRNAME}/*.txx ${SRC_DIRNAME}/CMakeLists.txt 
                  ${SRC_DIRNAME}/*doxygen.config)

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh 
                      ${SRC_DIRNAME}/ijkgenscalar.tar ${program_files} )

ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ../.. COMMAND zip
                      ${SRC_DIRNAME}/ijkgenscalar.zip ${program_files} )

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkgenscalar_doxygen.config)

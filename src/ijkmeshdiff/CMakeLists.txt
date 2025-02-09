CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(IJKMESHDIFF)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMESHDIFF_DIR "src/ijkmeshdiff")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${IJK_DIR}/lib")

ADD_COMPILE_OPTIONS(-Wall -Wno-unused-function)

SET(CMAKE_CXX_FLAGS "-std=c++20")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

ADD_EXECUTABLE(ijkmeshdiff ijkmeshdiff.cpp)

SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkmeshdiff DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkmeshdiff DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkmeshdiff DESTINATION "bin")
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKMESHDIFF_DIR}/ijkmcube.tar ${IJKMESHDIFF_DIR}/README ${IJKMESHDIFF_DIR}/INSTALL ${IJKMESHDIFF_DIR}/*.cpp ${IJKMESHDIFF_DIR}/*.h ${IJKMESHDIFF_DIR}/*.txx ${IJKMESHDIFF_DIR}/CMakeLists.txt )





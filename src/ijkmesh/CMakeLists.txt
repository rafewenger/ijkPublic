PROJECT(IJKMESH)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMESH_DIR "src/ijkmesh")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${IJK_DIR}/lib")

SET(CMAKE_CXX_FLAGS "-std=c++20")

ADD_EXECUTABLE(ijkmesh ijkmesh.cpp)

SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkmesh DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkmesh DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkmesh DESTINATION "bin")
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()

# OBSOLETE: SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
# OBSOLETE: INSTALL(TARGETS ijkmesh DESTINATION "bin")

SET(program_files ${IJKMESH_DIR}/README ${IJKMESH_DIR}/INSTALL 
                  ${IJKMESH_DIR}/*.cpp ${IJKMESH_DIR}/*.h 
                  ${IJKMESH_DIR}/*.tpp ${IJKMESH_DIR}/CMakeLists.txt 
                  ${IJKMESH_DIR}/*doxygen.config)

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh 
                      ${IJKMESH_DIR}/ijkmesh.tar ${program_files} )

ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ../.. COMMAND zip
                      ${IJKMESH_DIR}/ijkmesh.zip ${program_files} )

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkmesh_doxygen.config)

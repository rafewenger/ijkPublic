PROJECT(IJKMESHINFO)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED IJK_DIR)
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED IJK_DIR)

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMESHINFO_DIR "src/ijkmeshinfo")

#---------------------------------------------------------


INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${IJK_DIR}/lib")

SET(CMAKE_CXX_FLAGS "-std=c++11")

ADD_EXECUTABLE(ijkmeshinfo ijkmeshinfo.cpp ijkmeshinfo_compute.cpp
                           ijkmeshinfoIO.cpp)

SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkmeshinfo DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkmeshinfo DESTINATION "bin")
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkmeshinfo DESTINATION "bin")
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()

# OBSOLETE: SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
# OBSOLETE: INSTALL(TARGETS ijkmeshinfo DESTINATION "bin")

SET(program_files ${IJKMESHINFO_DIR}/README ${IJKMESHINFO_DIR}/INSTALL 
                  ${IJKMESHINFO_DIR}/*.cxx ${IJKMESHINFO_DIR}/*.h 
                  ${IJKMESHINFO_DIR}/*.txx ${IJKMESHINFO_DIR}/CMakeLists.txt 
                  ${IJKMESHINFO_DIR}/*doxygen.config)

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh 
                      ${IJKMESHINFO_DIR}/ijkmeshinfo.tar ${program_files} )
ADD_CUSTOM_TARGET(zip WORKING_DIRECTORY ../.. COMMAND zip
                      ${IJKMESHINFO_DIR}/ijkmeshinfo.zip ${program_files} )

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkmeshinfo_doxygen.config)

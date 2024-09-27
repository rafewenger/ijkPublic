PROJECT(IJKMCTABLE)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMCTABLE_DIR "src/ijkMCtable")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_LIBRARIES(expat)
ADD_EXECUTABLE(ijkgenMCtable ijkgenMCtable.cpp ijkgenMCpatch.cpp 
                             ijkMCtable.cpp ijkMCtable_poly.cpp 
                             ijkMCtable_xitIO.cpp clarkson_hull.c)
ADD_EXECUTABLE(ijkMCtable_invert ijkMCtable_invert.cpp
                                 ijkMCtable.cpp ijkMCtable_poly.cpp
                                 ijkMCtable_xitIO.cpp)
ADD_EXECUTABLE(ijkdiffMCtable ijkdiffMCtable.cpp ijkMCtable.cpp 
                              ijkMCtable_poly.cpp ijkMCtable_xitIO.cpp)
ADD_EXECUTABLE(ijkMCtableinfo ijkMCtableinfo.cpp ijkMCtable.cpp 
                              ijkMCtable_poly.cpp ijkMCtable_xitIO.cpp)
ADD_EXECUTABLE(ijktext2c ijktext2c.cpp)
ADD_EXECUTABLE(outambigtable outambigtable.cpp 
                             ijkMCtable_poly.cpp ijkMCtable_ambig.cpp)

SET(IJK_INSTALL_DIRA "${IJK_DIR}/bin")
SET(USER_SOFTWARE_DIR "/c/users/$ENV{USERNAME}/Software")
SET(IJK_INSTALL_DIRB "${USER_SOFTWARE_DIR}/bin")

IF (DEFINED IJK_INSTALL_DIR)
  SET(CMAKE_INSTALL_PREFIX ${IJK_INSTALL_DIR})
  INSTALL(TARGETS ijkgenMCtable DESTINATION "bin")
  INSTALL(TARGETS ijkMCtableinfo DESTINATION "bin")  
  INSTALL(TARGETS ijkdiffMCtable DESTINATION "bin")  
ELSEIF (EXISTS ${IJK_INSTALL_DIRA})
  SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
  INSTALL(TARGETS ijkgenMCtable DESTINATION "bin")
  INSTALL(TARGETS ijkMCtableinfo DESTINATION "bin")  
  INSTALL(TARGETS ijkdiffMCtable DESTINATION "bin")  
ELSEIF (EXISTS ${IJK_INSTALL_DIRB})
  SET(CMAKE_INSTALL_PREFIX ${USER_SOFTWARE_DIR})
  INSTALL(TARGETS ijkgenMCtable DESTINATION "bin")
  INSTALL(TARGETS ijkMCtableinfo DESTINATION "bin")  
  INSTALL(TARGETS ijkdiffMCtable DESTINATION "bin")  
ELSE()
  MESSAGE(WARNING " Warning: Unable to determine install directory.\n"
                  " Use -DIJK_INSTALL_DIR=... to specify install directory.\n"
                  " Create directory ${IJK_INSTALL_DIRA} and rerun cmake"
                  " if you want cmake to automatically use ${IJK_INSTALL_DIRA}"
                  " as the install directory.")
ENDIF()

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKMCTABLE_DIR}/ijktable.tar ${IJKMCTABLE_DIR}/README ${IJKMCTABLE_DIR}/INSTALL ${IJKMCTABLE_DIR}/*.cpp ${IJKMCTABLE_DIR}/*.h ${IJKMCTABLE_DIR}/*.tpp ${IJKMCTABLE_DIR}/*.c ${IJKMCTABLE_DIR}/CMakeLists.txt ${IJKMCTABLE_DIR}/ijktable.xsd ${IJKMCTABLE_DIR}/man/* )

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkMCtable_doxygen.config)

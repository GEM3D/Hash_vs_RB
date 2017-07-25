cmake_minimum_required(VERSION 2.8.11)

project(amrGem)

#disable in source build
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

message(STATUS "cmake is building on a ${CMAKE_SYSTEM_NAME} system")
set(dir ${PROJECT_SOURCE_DIR})
message("project source dir" ${dir})

set(CMAKE_FILES_DIRECTORY ${dir}/build)
message("Files Directory" ${CMAKE_FILES_DIRECTORY})

#define some colors (Unix-based systems)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()


# get the number of physical procesors 

message("project binary dir" ${EXECUTABLE_OUTPUT_PATH})
message("CMAKE SOURCE DIR:" ${CMAKE_SOURCE_DIR})
message("CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

# Out-of-Source build is required by amrGem, this if statement prevents
# in-source-builds 

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()


message("cmake Home Dir:" ${CMAKE_HOME_DIRECTORY})
message("cmake Build Dir:" ${CMAKE_BUILD_DIRECTORY})
message("cmake Binary Dir:" ${CMAKE_BINARY_DIR})
message("cmake source dir:" ${CMAKE_SOURCE_DIR})


############################################################################
#
#                    SECTION: FINDING-PACKAGE
#
#        1) HDF5 have native FindXXX.cmake's provided by cmake
#        2) For Zoltan and ParMETIS one needs to specify FindXXX.cmake 
#        3) These non-native to cmake modules are included in CMakeModules
#
############################################################################

#find HDF5
find_package(HDF5 ${HDF5_MIN_VERSION} REQUIRED)
if(HDF5_FOUND)
include_directories(${HDF5_INCLUDE_DIRS})
#message( "${Red}HDF5 Library Found: ${ColourReset}" ${HDF5_INCLUDE_DIRS})
message("${Red} HDF5 Library Found: ${ColourReset}")
endif()

###########################################################################
#
#                    SECTION: BUILD AND INSTALL
#
#        1) Build all the *.cpp and *.c files in "/src" directory
#        2) Link with external package libraries 
#        3) Install the final executable in "/bin"
#
############################################################################

include_directories(${dir}/src/include)
# dont forget to modify the address of the source files
file(GLOB OBJ "${dir}/src/*.cpp")
#set(CMAKE_BUILD_FILES_DIRECTORY ${dir})
#set(CMAKE_BUILD_DIRECTORY ${dir}/bin)      

set(EXECUTABLE_OUTPUT_PATH ${dir}/bin)

# Define all header files for linking, this is required such that one doesnt
# have to hardcode the directories for I/O operations and hence ensure
# portability

set(HEADER_FILES "${dir}/src/include/header.h" "${dir}/src/include/parse_stl.h")

message("header files" ${HEADER_FILES})

# insert -pg instead of -g for checking with gprof
# gprof bin/amrGem > out.text  
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++1y -pg -O3 -Wunused-variable ")
#    SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -Wall -std=gnu++0x" )
    SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -Wall -std=c++1y -Wunused-variable" )
endif()

# Generate Object Files
add_executable(amrGem ${OBJ} ${HEADER_FILES}) 

target_link_libraries(amrGem  ${HDF5_C_LIBRARIES})

#install(TARGETS amrGem DESTINATION ${dir/}bin)

############################################################################
#
#                    SECTION: ARCHIEVING
#
#	 1) it gets the date and appends that to the project name
#	 2) compresses the "/src" in a TGZ format
# 	 3) Puts the resulting file in the "/archieve" directory
#
############################################################################

MACRO (TODAY RESULT)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
         string(REGEX REPLACE "(..)/(..)/..(..).*" "\\2.\\1.\\3"
${RESULT} ${${RESULT}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\2.\\1.\\3"
${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

TODAY(DATE)
message(STATUS "Compilation date = ${DATE}")


# the following two lines specify the file type for binary and source
# generation separately

SET(CPACK_GENERATOR "TGZ")
SET(CPACK_SOURCE_GENERATOR "TGZ")

#INCLUDE(CPack)
# Do not rely on CPACK_xxx variables, they are not set till INCLUDE(CPack) is
# invoked, so there is no way to set them and as cpack to name the file as you
# like, I hard coded the project name and the date as a file_name to be used
# by cpack, this works fine 

#message("package name: " ${CPACK_PACKAGE_NAME})

set(CPACK_PACKAGE_FILE_NAME "amrGem-${DATE}")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "amrGem-${DATE}")

# set the directory you would like to be compressed, For under development
# projects it is normally the src directory that includes *.cpp and *.h files

SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${dir}/src;/")

SET(CPACK_PACKAGE_VENDOR "BSU")

# If you like to add any other files to archieve, need to use install file
# command, now in addition to source I am compresing the CMakeList.txt along
# with itimply installing it to the source directory. since src is the primary
# directory to archieve, make install command will copy this to src directory
# though, use the  following line if you would like to include
# CMakeList.txt in your archieve otherwise comment it out, one way to 
# get rid of this is to use an alias to remove that file from src after
# packing

install(FILES ${dir}/CMakeLists.txt DESTINATION ${dir}/src)

message("file name:" ${CPACK_PACKAGE_FILE_NAME})

#install(TARGETS  amrGem ARCHIVE  DESTINATION archieve)

# This variable (CPACK_OUTPUT_FILE_PREFIX)  will put the result of amrGem.date.tar.gz into the folder
# named archieve, this is only available in  in CMake versions >= 2.8.3 

set(CPACK_OUTPUT_FILE_PREFIX "${dir}/archieve")

#message("PACKAGE Install Dir:" ${CPACK_PACKAGE_INSTALL_DIRECTORY})

INCLUDE(CPack)

#
#  Doxygen Integrated to Generate Documentation for the Code
#
#

option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API
documentation" ON)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)

  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it
correctly")
  endif()

  add_custom_target (doc ALL
                     ${DOXYGEN_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
                     )

endif()




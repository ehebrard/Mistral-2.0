cmake_minimum_required(VERSION 3.10)
project(schedcl)
include(FetchContent)

# Parser lib
FetchContent_Declare(xcsp3 GIT_REPOSITORY https://github.com/xcsp3team/XCSP3-CPP-Parser)
FetchContent_MakeAvailable(xcsp3)

set(XCSP3_INCLUDE ${xcsp3_SOURCE_DIR}/include)
include_directories(${XCSP3_INCLUDE})
file(GLOB_RECURSE XCSP3_SOURCES ${xcsp3_SOURCE_DIR}/src/*.cc)

find_package(LibXml2 REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-mtune=native -march=native")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
if (BUILD_TYPE STREQUAL "debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
    message("Building for debug")
elseif (BUILD_TYPE STREQUAL "release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -flto")
    message("Building for release")
else ()
    message(SEND_ERROR "Unknown build type ${CMAKE_BUILD_TYPE}")
endif ()

file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/lib/*.cpp)
include_directories(${CMAKE_SOURCE_DIR}/src/include/)
include_directories(${CMAKE_SOURCE_DIR}/tools/tclap/include/)

file(GLOB_RECURSE EXECUTABLES ${CMAKE_SOURCE_DIR}/examples/src/*.cpp)
message("Generating targets for the following cpp files:")
foreach(EXEC ${EXECUTABLES})
    get_filename_component(NAME ${EXEC} NAME_WLE)
    message(\t${EXEC}\ ->\ target:\ ${NAME})
    add_executable(${NAME} ${EXEC} ${SOURCES} ${XCSP3_SOURCES})
    target_link_libraries(${NAME} PRIVATE ${LIBXML2_LIBRARIES})
endforeach()

cmake_minimum_required (VERSION 2.6)

project (rbdli)
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++11 -O3 -DEIGEN_NO_DEBUG -w")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -DEIGEN_NO_DEBUG -w")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${rbdli_SOURCE_DIR})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${rbdli_SOURCE_DIR}/../)
SET(BOOST_MIN_VERSION "1.54.0")
SET(PYTHON_MIN_VERSION "2.7.0")
include(FindPkgConfig)

############################# LOAD RBDL #############################

if(PKG_CONFIG_FOUND)
    pkg_check_modules(RBDL rbdl)
    if (NOT RBDL_FOUND)
        message(FATAL_ERROR "RBDL not found")
    endif()
    include_directories(${RBDL_INCLUDE_DIRS})
	link_directories(${RBDL_LIBRARY_DIRS}) 
	message("RBDL inc dirs ${RBDL_INCLUDE_DIRS}") 
	message("RBDL lib dirs ${RBDL_LIBRARY_DIRS}")  
endif()

############################# LOAD URDFREADER #############################
if(PKG_CONFIG_FOUND)
    pkg_check_modules(URDFREADER urdfreader)
    if(NOT URDFREADER_FOUND)
       message(FATAL_ERROR "URDFREADER could not be found")
    endif()
    include_directories(${URDFREADER_INCLUDE_DIRS})
    link_directories(${URDFREADER_LIBRARY_DIRS})    
endif()

############################# LOAD BOOST #############################
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS python system thread timer)
if (NOT Boost_FOUND)
    set(Boost_INCLUDE_DIRS "${BOOST_INCLUDE_DIRS}")
    IF ( NOT Boost_INCLUDE_DIRS )
        message(FATAL_ERROR "Please point the environment variable BOOST_INCLUDE_DIRS to the include directory of Boost")
    ENDIF()
endif ()
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

############################# LOAD EIGEN #############################
if(PKG_CONFIG_FOUND)
    pkg_check_modules(EIGEN eigen3)
    if(NOT EIGEN_FOUND)
       message(FATAL_ERROR "EIGEN could not be found")
    endif()
    include_directories(${EIGEN_INCLUDE_DIRS})
    link_directories(${EIGEN_LIBRARY_DIRS})    
endif()

add_library(rbdl_interface SHARED 
            rbdl_interface.cpp)
            
  
            
target_link_libraries (rbdl_interface
                       ${RBDL_LIBRARIES}                       
                       ${EIGEN_LIBRARIES}
                       ${Boost_LIBRARIES}
                       ${URDFREADER_LIBRARIES})
                       
set(pkg_conf_file_in "${CMAKE_CURRENT_SOURCE_DIR}/rbdl_interface.pc.in")
set(pkg_conf_file_out "${CMAKE_CURRENT_BINARY_DIR}/rbdl_interface.pc")
configure_file("${pkg_conf_file_in}" "${pkg_conf_file_out}" @ONLY)

install(TARGETS rbdl_interface DESTINATION lib)
install(FILES rbdl_interface.hpp DESTINATION include/rbdl_interface)
install(FILES "${pkg_conf_file_out}" DESTINATION lib/pkgconfig/ COMPONENT pkgconfig)

#pkg_version=1.4.1
cmake_minimum_required(VERSION 2.8.3)
project(momo)

###################
## find packages ##
###################
find_package(mrt_cmake_modules REQUIRED)
include(UseMrtStdCompilerFlags)
include(UseMrtAutoTarget)

include(GatherDeps)
#remove libs, which cannot be found automatically
#list(REMOVE_ITEM DEPENDEND_PACKAGES <package name 1> <package name 2> ...)
find_package(AutoDeps REQUIRED COMPONENTS ${DEPENDEND_PACKAGES})

#manually resolve removed dependend packages
#find_package(...)
############################
## read source code files ##
############################
file(GLOB_RECURSE PROJECT_HEADER_FILES_INC RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "include/*.h" "include/*.hpp")
file(GLOB PROJECT_SOURCE_FILES_INC RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "src/*.h" "src/*.hpp")
file(GLOB PROJECT_SOURCE_FILES_SRC RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "src/*.cpp")
if (PROJECT_SOURCE_FILES_SRC)
	set(LIBRARY_NAME ${PROJECT_NAME})
endif()

###################################
## catkin specific configuration ##
###################################
# The catkin_package macro generates cmake config files for your package
# Declare things to be passed to dependent projects
# INCLUDE_DIRS: uncomment this if you package contains header files
# LIBRARIES: libraries you create in this project that dependent projects also need
# CATKIN_DEPENDS: catkin_packages dependent projects also need
# DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS include ${mrt_EXPORT_INCLUDE_DIRS}
  LIBRARIES ${LIBRARY_NAME} ${mrt_EXPORT_LIBRARIES}
  CATKIN_DEPENDS ${catkin_EXPORT_DEPENDS}
)

###########
## Build ##
###########
# Add include and library directories
include_directories(
  include/${LIBRARY_NAME}
  ${mrt_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
)

link_directories(
  ${mrt_LIBRARY_DIRS}
)

if (PROJECT_SOURCE_FILES_SRC)
	# Declare a cpp library
	add_library(${LIBRARY_NAME}
	  ${PROJECT_HEADER_FILES_INC}
	  ${PROJECT_SOURCE_FILES_INC}
	  ${PROJECT_SOURCE_FILES_SRC}
	)
	# Specify libraries to link a library or executable target against
	target_link_libraries(${LIBRARY_NAME}
	  ${catkin_LIBRARIES}
	  ${mrt_LIBRARIES}
	)
endif()

#############
## Install ##
#############
if (TARGET ${PROJECT_NAME})
	# Mark library for installation
	install(TARGETS ${PROJECT_NAME}
	  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
	  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
	  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
	)
endif()
# Mark c++ header files for installation
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
#install(FILES
#  res/test.png
#  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
#)

#############
## Testing ##
#############
# Add gtest based cpp test target and link libraries
if (CATKIN_ENABLE_TESTING)
	file(GLOB PROJECT_TEST_FILES_SRC RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "test/*.cpp")
	foreach(PROJECT_TEST_FILE_SRC ${PROJECT_TEST_FILES_SRC})
		get_filename_component(PROJECT_TEST_NAME ${PROJECT_TEST_FILE_SRC} NAME_WE)
		catkin_add_gtest(${PROJECT_TEST_NAME}-test ${PROJECT_TEST_FILE_SRC})
		target_link_libraries(${PROJECT_TEST_NAME}-test ${LIBRARY_NAME} ${catkin_LIBRARIES} ${mrt_LIBRARIES} gtest_main)
	endforeach()
endif()

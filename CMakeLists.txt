cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(Template
  VERSION 0.1
  LANGUAGES CXX C
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")

include(FetchContent)
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

# list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(FeatureSummary)

# include(CodeFormat)

find_package(ProXPDE REQUIRED)

add_executable(main)
target_sources(main
	PRIVATE main.cpp
	# PUBLIC ${HEADER_FILES}
)
target_link_libraries(main
  PUBLIC
    ProXPDE::proxpde
    fmt::fmt
)

add_executable(mysolver)
target_sources(mysolver
	PRIVATE mysolver.cpp
	# PUBLIC ${HEADER_FILES}
)
target_link_libraries(mysolver
  PUBLIC
    ProXPDE::proxpde
    fmt::fmt
)

feature_summary(WHAT ALL)

set(SRC_FILES
  main.cpp
  mysolver.cpp
)
# set(HEADER_FILES
#   header.hpp
# )

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  add_custom_target(format
    COMMAND ${CLANGFORMAT}
    --verbose -i
    ${SRC_FILES} ${HEADER_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
else()
  add_custom_target(format
    COMMAND echo "clang-format could not be found"
  )
endif()


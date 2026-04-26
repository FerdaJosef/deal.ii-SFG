##
#  CMake script for the step-3 tutorial program:
##

set(MODEL_DIR "${CMAKE_SOURCE_DIR}/source/Model3D_2vars")

file(COPY ${MODEL_DIR}/double_ditch.prm
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${MODEL_DIR}
    ${MODEL_DIR}/AceGen
    ${CMAKE_SOURCE_DIR}/include
)

# Set the name of the project and target:
set(TARGET "main")

#file(GLOB_RECURSE TARGET_SRC  "source/Model2D_3vars/*.cc")

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#    file(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#    file(GLOB_RECURSE TARGET_INC  "include/*.h")
#    set(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
set(TARGET_SRC
    ${MODEL_DIR}/model.cc
    ${MODEL_DIR}/InitialValues.cc
    ${MODEL_DIR}/RandomField.cc
    main.cc
    )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.7.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

deal_ii_initialize_cached_variables()
project(${TARGET})
deal_ii_invoke_autopilot()


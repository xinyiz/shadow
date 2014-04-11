
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(voxelizer)

file(GLOB_RECURSE HEADER_CODE ${voxelizer_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE SRC_CODE ${voxelizer_SOURCE_DIR}/source/*.cpp)

ADD_EXECUTABLE(voxelizer ${SRC_CODE} ${HEADER_CODE})
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIRS})
target_link_libraries(voxelizer ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${voxelizer_SOURCE_DIR}/lib/libvecmath.a)

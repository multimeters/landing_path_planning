cmake_minimum_required(VERSION 3.1)
project(VFRRT)
set(CMAKE_BUILD_TYPE Debug)  # 确保是 Debug 模式
# 设置Eigen的路径
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# 手动设置OMPL的路径
set(OMPL_INCLUDE_DIRS /usr/local/include/ompl-1.6)
set(OMPL_LIBRARIES /usr/local/lib/libompl.so) # 或者你系统中实际的OMPL库文件

include_directories(${OMPL_INCLUDE_DIRS})

add_executable(vf_rrt main_vf_down_sample16.cpp)

# 链接OMPL库
target_link_libraries(vf_rrt ${OMPL_LIBRARIES})

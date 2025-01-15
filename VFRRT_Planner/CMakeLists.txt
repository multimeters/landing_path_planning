cmake_minimum_required(VERSION 3.1)
project(VFRRTPlanner)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)

# 手动指定 OMPL 的安装路径
# 请根据您的实际安装路径进行修改
set(OMPL_INCLUDE_DIR "/usr/local/include/ompl-1.6")
set(OMPL_LIB_DIR "/usr/local/lib")

# 包含 OMPL 的头文件
include_directories(${OMPL_INCLUDE_DIR})

# 添加链接目录
link_directories(${OMPL_LIB_DIR})

# 查找 OMPL 的库文件
find_library(OMPL_LIB ompl PATHS ${OMPL_LIB_DIR} REQUIRED)

# 检查是否找到 OMPL 库
if(NOT OMPL_LIB)
    message(FATAL_ERROR "OMPL library not found in ${OMPL_LIB_DIR}")
endif()

# 查找 Eigen 库
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 包含 Eigen 的头文件
include_directories(${EIGEN3_INCLUDE_DIR})

# 如果使用 `find_package(Eigen3 REQUIRED)`，通常会自动设置 `EIGEN3_INCLUDE_DIR`

# 添加可执行文件
add_executable(VFRRTPlanner src/VFRRTPlanner_multi_3.cpp)

# 链接 OMPL 库
target_link_libraries(VFRRTPlanner ${OMPL_LIB} Eigen3::Eigen)
# 如果 OMPL 依赖于其他库（如 Boost），请确保链接这些库
# 例如：
# find_package(Boost REQUIRED COMPONENTS system filesystem)
# include_directories(${Boost_INCLUDE_DIRS})
# target_link_libraries(VFRRTPlanner ${OMPL_LIB} ${Boost_LIBRARIES} Eigen3::Eigen)

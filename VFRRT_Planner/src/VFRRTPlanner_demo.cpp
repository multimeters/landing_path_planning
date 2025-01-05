#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/config.h>
#include <fstream>
#include <iostream>

namespace ob = ompl::base;
namespace og = ompl::geometric;

// 定义环境中的障碍物（简单示例）
bool isStateValid(const ob::State *state)
{
    const auto *pos = state->as<ob::RealVectorStateSpace::StateType>();
    double x = pos->values[0];
    double y = pos->values[1];

    // 简单的圆形障碍物
    double obstacle_center_x = 5.0;
    double obstacle_center_y = 5.0;
    double obstacle_radius = 1.0;

    double dist = std::sqrt((x - obstacle_center_x)*(x - obstacle_center_x) + (y - obstacle_center_y)*(y - obstacle_center_y));
    return dist > obstacle_radius;
}

int main()
{
    // 1. 定义状态空间 (2D)
    auto space(std::make_shared<ob::RealVectorStateSpace>(2));

    // 2. 设置空间的边界
    ob::RealVectorBounds bounds(2);
    bounds.setLow(-10);
    bounds.setHigh(10);
    space->setBounds(bounds);

    // 3. 创建SimpleSetup
    og::SimpleSetup ss(space);

    // 4. 设置状态有效性检查
    ss.setStateValidityChecker(isStateValid);

    // 5. 设置起始和目标状态
    ob::ScopedState<> start(space);
    start->as<ob::RealVectorStateSpace::StateType>()->values[0] = -8.0;
    start->as<ob::RealVectorStateSpace::StateType>()->values[1] = -8.0;

    ob::ScopedState<> goal(space);
    goal->as<ob::RealVectorStateSpace::StateType>()->values[0] = 8.0;
    goal->as<ob::RealVectorStateSpace::StateType>()->values[1] = 8.0;

    ss.setStartAndGoalStates(start, goal);

    // 6. 使用RRT规划器
    ob::PlannerPtr planner(new og::RRT(ss.getSpaceInformation()));
    ss.setPlanner(planner);

    // 7. 进行规划
    ob::PlannerStatus solved = ss.solve(5.0);

    if (solved)
    {
        std::cout << "找到路径!" << std::endl;

        // 8. 获取路径
        og::PathGeometric path = ss.getSolutionPath();
        path.interpolate(); // 插值以获得更平滑的路径

        // 9. 保存路径到文件
        std::ofstream path_file("data/path.txt");
        for (std::size_t i = 0; i < path.getStateCount(); ++i)
        {
            const auto *state = path.getState(i)->as<ob::RealVectorStateSpace::StateType>();
            path_file << state->values[0] << " " << state->values[1] << "\n";
        }
        path_file.close();

        // 10. 生成向量场并保存到文件（简单示例：梯度场）
        std::ofstream vf_file("data/vector_field.txt");
        double resolution = 1.0;
        for(double x = -10; x <= 10; x += resolution)
        {
            for(double y = -10; y <= 10; y += resolution)
            {
                // 简单的向量场示例：指向目标的方向
                double vx = 8.0 - x;
                double vy = 8.0 - y;
                double norm = std::sqrt(vx*vx + vy*vy);
                if(norm != 0)
                {
                    vx /= norm;
                    vy /= norm;
                }
                vf_file << x << " " << y << " " << vx << " " << vy << "\n";
            }
        }
        vf_file.close();

        // 可视化路径
        ss.simplifySolution();
        ss.getSolutionPath().printAsMatrix(std::cout);
    }
    else
    {
        std::cout << "未找到路径。" << std::endl;
    }

    return 0;
}

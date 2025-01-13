#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/VFRRT.h>
#include <ompl/util/RandomNumbers.h>
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <ompl/base/objectives/VFUpstreamCriterionOptimizationObjective.h>
#include <filesystem> // 添加此头文件以使用 filesystem

// 使用OMPL的命名空间简化代码
namespace ob = ompl::base;
namespace og = ompl::geometric;
namespace fs = std::filesystem;

// 定义向量场类型，使用std::function封装接受OMPL状态并返回Eigen向量的函数
using VectorField = std::function<Eigen::VectorXd(const ob::State *)>;

/**
 * @brief 从文件中加载矩阵数据
 * 
 * @param filename 文件路径
 * @param rows 矩阵行数
 * @param cols 矩阵列数
 * @return Eigen::MatrixXd 加载后的矩阵
 */
Eigen::MatrixXd loadMatrixFromFile(const std::string &filename, int rows, int cols)
{
    Eigen::MatrixXd matrix(rows, cols);
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return Eigen::MatrixXd(); // 返回空矩阵
    }

    // 按行逐元素读取数据
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (!(file >> matrix(i, j)))
            {
                std::cerr << "读取文件时遇到错误: " << filename << std::endl;
                return Eigen::MatrixXd(); // 返回空矩阵
            }
        }
    }

    file.close();
    return matrix;
}

/**
 * @brief 将路径保存到文件中
 * 
 * @param path 规划得到的路径
 * @param filename 保存的文件名
 */
void savePathToFile(const og::PathGeometric &path, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件进行写入: " << filename << std::endl;
        return;
    }

    // 遍历路径中的每一个状态并写入文件
    for (size_t i = 0; i < path.getStateCount(); ++i)
    {
        const auto *state = path.getState(i)->as<ob::RealVectorStateSpace::StateType>();
        file << state->values[0] << " " << state->values[1] << std::endl;
    }

    file.close();
}

/**
 * @brief 从文件中读取采样点
 * 
 * @param filename 文件路径
 * @param points 存储读取到的点的向量
 * @return true 读取成功
 * @return false 读取失败
 */
bool loadSamplePoints(const std::string &filename, std::vector<std::pair<double, double>> &points)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开采样点文件: " << filename << std::endl;
        return false;
    }

    double x, y;
    while (file >> x >> y)
    {
        points.emplace_back(x, y);
    }

    file.close();
    return true;
}

/**
 * @brief 将向量场保存到文件中，便于可视化或调试
 * 
 * @param vectorField 向量场函数
 * @param filename 保存的文件名
 */
void saveVectorFieldToFile(const VectorField &vectorField, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件进行写入: " << filename << std::endl;
        return;
    }

    // 定义采样点数
    const int numPointsX = 512;  // X轴采样点数
    const int numPointsY = 256;  // Y轴采样点数（根据向量场矩阵大小调整）

    // 计算每个采样点的步长
    const double stepX = 1.0;    // X轴步长为1米
    const double stepY = 2.0;    // Y轴步长为2米

    // 创建一个二维状态空间用于采样
    auto space = std::make_shared<ob::RealVectorStateSpace>(2);
    ob::ScopedState<ob::RealVectorStateSpace> state(space);

    // 遍历所有采样点，计算并保存向量场的值
    for (int i = 0; i < numPointsY; ++i)
    {
        double y = i * stepY;
        for (int j = 0; j < numPointsX; ++j)
        {
            double x = j * stepX;

            state->values[0] = x;  // X坐标
            state->values[1] = y;  // Y坐标

            Eigen::VectorXd field = vectorField(state.get());

            file << x << " " << y << " " << field[0] << " " << field[1] << std::endl;
        }
    }

    file.close();
}

int main()
{
    // 定义二维实数向量状态空间
    auto space = std::make_shared<ob::RealVectorStateSpace>(2);

    // 设置状态空间的边界
    ob::RealVectorBounds bounds(2);
    bounds.setLow(0, 0.0);     // X轴下界为0米
    bounds.setHigh(0, 512.0);  // X轴上界为512米
    bounds.setLow(1, 0.0);     // Y轴下界为0米
    bounds.setHigh(1, 512.0);  // Y轴上界为500米
    space->as<ob::RealVectorStateSpace>()->setBounds(bounds);

    // 创建空间信息对象
    auto si = std::make_shared<ob::SpaceInformation>(space);

    // 从文件中加载向量场数据
    Eigen::MatrixXd u_combined = loadMatrixFromFile("/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/u_vf_normalized.txt", 256, 512);
    Eigen::MatrixXd v_combined = loadMatrixFromFile("/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/v_vf_normalized.txt", 256, 512);

    // 检查向量场数据是否加载成功
    if (u_combined.size() == 0 || v_combined.size() == 0)
    {
        std::cerr << "向量场数据加载失败。" << std::endl;
        return -1;
    }

    // 定义向量场函数，根据状态返回对应的向量值
    VectorField vectorField = [&u_combined, &v_combined](const ob::State *state) -> Eigen::VectorXd {
        const auto *realState = state->as<ob::RealVectorStateSpace::StateType>();
        double x = realState->values[0] / 1.0;  // 将米转换为下采样后的矩阵索引
        double y = realState->values[1] / 2.0;  // 将米转换为下采样后的矩阵索引

        // 将坐标映射到矩阵索引
        int i = static_cast<int>(std::round(y)); // 行索引，对应Y轴
        int j = static_cast<int>(std::round(x)); // 列索引，对应X轴

        // 确保索引在有效范围内
        i = std::min(std::max(i, 0), static_cast<int>(u_combined.rows() - 1));
        j = std::min(std::max(j, 0), static_cast<int>(u_combined.cols() - 1));

        Eigen::VectorXd field(2);
        field[0] = u_combined(i, j);
        field[1] = v_combined(i, j);
        return field;
    };

    // 可选：将向量场保存到文件，便于可视化或调试
    saveVectorFieldToFile(vectorField, "vector_Field.txt");

    // 初始化SimpleSetup对象，管理规划过程
    og::SimpleSetup ss(space);

    // 设置状态有效性检查器，假设所有状态都是有效的
    ss.setStateValidityChecker(std::make_shared<ob::AllValidStateValidityChecker>(si));

    // 定义起始状态
    ob::ScopedState<> start(space);
    start[0] = 400.0;  // 起始点的X坐标为400米
    start[1] = 350.0;  // 起始点的Y坐标为350米

    // 设置起始状态
    ss.setStartState(start);

    // 加载采样点作为目标点
    std::vector<std::pair<double, double>> samplePoints;
    std::string samplePointsFile = "/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/sampled_coastline.txt"; // 请确保路径正确
    if (!loadSamplePoints(samplePointsFile, samplePoints))
    {
        std::cerr << "加载采样点失败。" << std::endl;
        return -1;
    }

    if (samplePoints.empty())
    {
        std::cerr << "没有读取到任何采样点。" << std::endl;
        return -1;
    }

    // 确保 "path" 文件夹存在
    std::string pathDir = "path";
    if (!fs::exists(pathDir))
    {
        if (!fs::create_directory(pathDir))
        {
            std::cerr << "无法创建目录: " << pathDir << std::endl;
            return -1;
        }
    }

    // 设置规划器的参数
    double exploration = 0.1;                  // 探索参数
    double initial_lambda = 100.1;             // 初始lambda参数
    unsigned int initial_lambda_samples = 1000000; // 初始lambda采样数

    // 在循环外创建并设置规划器
    auto planner = std::make_shared<og::VFRRT>(ss.getSpaceInformation(), vectorField, exploration, initial_lambda, initial_lambda_samples);
    ss.setPlanner(planner);
    ss.setup(); // 初始设置

    // 遍历每个采样点，设置为目标并进行路径规划
    int pathIndex = 0;
    for (const auto &point : samplePoints)
    {
        double goalX = point.first;
        double goalY = point.second;

        // 定义目标状态
        ob::ScopedState<> goal(space);
        goal[0] = goalX;     // 目标点的X坐标
        goal[1] = goalY;     // 目标点的Y坐标

        std::cout << "设置新的目标为: (" << goalX << ", " << goalY << ")\n"; // 调试输出

        // 设置新的目标状态
        ss.setGoalState(goal);
        ss.setup(); // 更新规划器

        // 清除之前的解决方案
        ss.getProblemDefinition()->clearSolutionPaths();

        // 尝试在指定时间内解决规划问题（例如，400秒）
        ob::PlannerStatus solved = ss.solve(400.0);

        if (solved)
        {
            if (solved == ob::PlannerStatus::EXACT_SOLUTION)
                std::cout << "路径 " << pathIndex << ": 找到精确解。\n";
            else
                std::cout << "路径 " << pathIndex << ": 找到近似解。\n";

            // 获取并插值解决方案路径
            og::PathGeometric path = ss.getSolutionPath();
            path.interpolate();

            // 使用向量场计算总上游成本
            auto upstream = std::make_shared<ob::VFUpstreamCriterionOptimizationObjective>(ss.getSpaceInformation(), vectorField);
            double totalUpstreamCost = path.cost(upstream).value();
            std::cout << "路径 " << pathIndex << ": 总上游成本: " << totalUpstreamCost << "\n";

            // 生成文件名，例如 "path_0.txt", "path_1.txt", ...
            std::string filename = pathDir + "/path_" + std::to_string(pathIndex) + ".txt";
            std::cout << "路径 " << pathIndex << " 将保存到: " << filename << std::endl;
            savePathToFile(path, filename);
        }
        else
        {
            std::cout << "路径 " << pathIndex << ": 未找到解决方案。\n";
        }

        pathIndex++;
    }

    std::cout << "所有路径规划完成。" << std::endl;
    return 0;
}

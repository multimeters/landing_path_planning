#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/VFRRT.h>
#include <ompl/util/RandomNumbers.h>
#include <ompl/base/objectives/VFUpstreamCriterionOptimizationObjective.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/tools/config/MagicConstants.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <filesystem>
#include <thread>
#include <future>
#include <mutex>
#include <chrono>  // 用于测量程序运行时间
#include <getopt.h>  // 用于解析命令行参数
// 使用 OMPL 的命名空间简化代码
namespace ob = ompl::base;
namespace og = ompl::geometric;
namespace fs = std::filesystem;

// 定义向量场类型，使用 std::function 封装
using VectorField = std::function<Eigen::VectorXd(const ob::State *)>;

/**
 * @brief 从文件中加载矩阵数据
 *
 * @param filename 文件路径
 * @param rows 矩阵行数
 * @param cols 矩阵列数
 * @return Eigen::MatrixXd 加载后的矩阵；若出错则返回空矩阵
 */
Eigen::MatrixXd loadMatrixFromFile(const std::string &filename, int rows, int cols)
{
    Eigen::MatrixXd matrix(rows, cols);
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "[Error] 无法打开文件: " << filename << std::endl;
        return Eigen::MatrixXd(); // 返回空矩阵
    }

    // 按行逐元素读取数据
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (!(file >> matrix(i, j)))
            {
                std::cerr << "[Error] 读取文件时遇到错误: " << filename << std::endl;
                return Eigen::MatrixXd(); // 返回空矩阵
            }
        }
    }
    file.close();

    // 若一切正常，返回 matrix
    return matrix;
}

/**
 * @brief 将路径保存到文件
 *
 * @param path 规划得到的路径
 * @param filename 保存的文件名
 */
void savePathToFile(const og::PathGeometric &path, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "[Error] 无法打开文件进行写入: " << filename << std::endl;
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
 * @brief 从文件中读取采样点 (x, y)
 *
 * @param filename 文件路径
 * @param points 存储读取到的点 (x,y) 的向量
 * @return true 如果读取成功
 * @return false 如果失败
 */
bool loadSamplePoints(const std::string &filename, std::vector<std::pair<double, double>> &points)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "[Error] 无法打开采样点文件: " << filename << std::endl;
        return false;
    }

    double x, y;
    while (file >> x >> y)
    {
        points.emplace_back(x, y);
    }
    file.close();

    return true; // 读取成功后记得返回 true
}

/**
 * @brief 将向量场保存到文件，便于可视化或调试
 *
 * @param vectorField 向量场函数
 * @param filename 保存的文件名
 */
void saveVectorFieldToFile(const VectorField &vectorField, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "[Error] 无法打开文件进行写入: " << filename << std::endl;
        return;
    }

    // 定义采样点数
    const int numPointsX = 512;  // X轴采样点数
    const int numPointsY = 256;  // Y轴采样点数（根据向量场矩阵大小调整）

    // 计算每个采样点的步长
    const double stepX = 1.0;    // X轴步长为 1 米
    const double stepY = 2.0;    // Y轴步长为 2 米

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
            state->values[0] = x;
            state->values[1] = y;

            Eigen::VectorXd field = vectorField(state.get());
            file << x << " " << y << " " << field[0] << " " << field[1] << std::endl;
        }
    }
    file.close();
}

/**
 * @brief 用于在多线程之间传递规划信息的结构体
 */
struct PathPlanningResult
{
    bool foundSolution{false};   // 是否找到解
    double upstreamCost{0.0};    // 上游成本
    int pathIndex{-1};           // 路径索引

    // ======= 修改部分：增加一个标识，表示是否是 "EXACT_SOLUTION" =======
    bool isExactSolution{false}; // 是否是精确解
};

/**
 * @brief 规划并保存路径的函数 (在线程中执行)
 *
 * @param vectorField 向量场函数
 * @param point 目标采样点 (goalX, goalY)
 * @param pathDir 保存路径的目录
 * @param index 路径的索引
 * @param exploration VFRRT 的探索参数
 * @param initial_lambda 初始 lambda 参数
 * @param initial_lambda_samples 初始 lambda 采样数
 * @param u_combined 预加载的 u 向量场
 * @param v_combined 预加载的 v 向量场
 * @param coutMutex 互斥锁用于同步输出
 * @param logFile 用于将日志信息写入文件
 * @return PathPlanningResult 包含是否找到解、上游成本、路径索引、是否为精确解
 */
PathPlanningResult planAndSavePath(
    const VectorField &vectorField,
    const std::pair<double, double> &point,
    const std::string &pathDir,
    int index,
    double exploration,
    double initial_lambda,
    unsigned int initial_lambda_samples,
    const Eigen::MatrixXd &u_combined,
    const Eigen::MatrixXd &v_combined,
    std::mutex &coutMutex,
    std::ofstream &logFile)
{
    PathPlanningResult result;
    result.pathIndex = index;

    double goalX = point.first;
    double goalY = point.second;

    {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "路径 " << index << " 开始规划。\n"
                  << "目标: (" << goalX << ", " << goalY << ")\n";
        logFile << "[Info] 路径 " << index << " 开始规划, 目标:("
                << goalX << ", " << goalY << ")\n";
    }

    // ===== 1) 创建状态空间并设置边界 =====
    auto space = std::make_shared<ob::RealVectorStateSpace>(2);

    ob::RealVectorBounds bounds(2);
    bounds.setLow(0, 0.0);
    bounds.setHigh(0, 512.0);
    bounds.setLow(1, 0.0);
    bounds.setHigh(1, 400.0);
    space->as<ob::RealVectorStateSpace>()->setBounds(bounds);

    auto si = std::make_shared<ob::SpaceInformation>(space);

    // ===== 2) 初始化 SimpleSetup =====
    og::SimpleSetup ss(space);
    ss.setStateValidityChecker(std::make_shared<ob::AllValidStateValidityChecker>(si));

    // ===== 3) 设置起点和目标点 =====
    ob::ScopedState<> start(space);
    start[0] = 400.0;  // 可根据需要修改
    start[1] = 350.0;
    ss.setStartState(start);

    ob::ScopedState<> goal(space);
    goal[0] = goalX;
    goal[1] = goalY;
    ss.setGoalState(goal);

    // ===== 4) 创建并设置 VFRRT 规划器 =====
    auto planner = std::make_shared<og::VFRRT>(
        ss.getSpaceInformation(),
        vectorField,
        exploration,
        initial_lambda,
        initial_lambda_samples);

    ss.setPlanner(planner);
    ss.setup();

    // ===== 5) 尝试在 400 秒内规划 =====
    ob::PlannerStatus solved = ss.solve(500.0);

    if (solved)
    {
        // 找到解
        result.foundSolution = true;

        // ======= 修改部分：判断是否是 EXACT_SOLUTION =======
        if (solved == ob::PlannerStatus::EXACT_SOLUTION)
        {
            result.isExactSolution = true;
        }

        // 获取并插值路径
        og::PathGeometric path = ss.getSolutionPath();
        path.interpolate();

        // 计算上游成本
        auto upstream = std::make_shared<ob::VFUpstreamCriterionOptimizationObjective>(
            ss.getSpaceInformation(), vectorField);
        double totalUpstreamCost = path.cost(upstream).value();
        result.upstreamCost = totalUpstreamCost;

        {
            std::lock_guard<std::mutex> lock(coutMutex);

            if (result.isExactSolution)
            {
                std::cout << "路径 " << index << ": 找到精确解。\n";
                logFile << "[Info] 路径 " << index << " 找到精确解。\n";
            }
            else
            {
                std::cout << "路径 " << index << ": 找到近似解。\n";
                logFile << "[Info] 路径 " << index << " 找到近似解。\n";
            }

            std::cout << "路径 " << index << " 上游成本: " << totalUpstreamCost << "\n";
            logFile << "[Info] 路径 " << index << " 上游成本: " << totalUpstreamCost << "\n";

            // 保存到文本文件
            std::string filename = pathDir + "/path_" + std::to_string(index) + ".txt";
            std::cout << "路径 " << index << " 将保存到: " << filename << std::endl;
            logFile << "[Info] 路径 " << index << " 将保存到: " << filename << std::endl;
            savePathToFile(path, filename);
        }
    }
    else
    {
        // 没找到解
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "路径 " << index << ": 未找到解决方案。\n";
        logFile << "[Warning] 路径 " << index << ": 未找到解决方案。\n";
    }

    // 标记规划结束
    {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "路径 " << index << " 规划完成。\n";
        logFile << "[Info] 路径 " << index << " 规划完成。\n";
    }

    return result;
}

int main()
{
    // 记录程序开始时间
    auto programStart = std::chrono::steady_clock::now();

    // ===== 1) 读取向量场数据 =====
    Eigen::MatrixXd u_combined = loadMatrixFromFile(
        "/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/u_combined.txt",
        256, 512);
    Eigen::MatrixXd v_combined = loadMatrixFromFile(
        "/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/v_combined.txt",
        256, 512);

    if (u_combined.size() == 0 || v_combined.size() == 0)
    {
        std::cerr << "[Error] 向量场数据加载失败。" << std::endl;
        return -1;
    }

    // ===== 2) 定义向量场函数 =====
    VectorField vectorField = [&u_combined, &v_combined](const ob::State *state) -> Eigen::VectorXd {
        const auto *realState = state->as<ob::RealVectorStateSpace::StateType>();
        double x = realState->values[0] / 1.0;  // 将米转换为列索引
        double y = realState->values[1] / 2.0;  // 将米转换为行索引

        int i = static_cast<int>(std::round(y)); // 行
        int j = static_cast<int>(std::round(x)); // 列

        // 确保索引不越界
        i = std::min(std::max(i, 0), (int)u_combined.rows() - 1);
        j = std::min(std::max(j, 0), (int)u_combined.cols() - 1);

        Eigen::VectorXd field(2);
        field[0] = u_combined(i, j);
        field[1] = v_combined(i, j);
        return field;
    };

    // 可选：将向量场保存到文件，便于可视化或调试
    saveVectorFieldToFile(vectorField, "vector_Field.txt");

    // ===== 3) 读取采样点 =====
    std::vector<std::pair<double, double>> samplePoints;
    std::string samplePointsFile = "/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/sampled_coastline.txt";
    if (!loadSamplePoints(samplePointsFile, samplePoints))
    {
        std::cerr << "[Error] 加载采样点失败。" << std::endl;
        return -1;
    }
    if (samplePoints.empty())
    {
        std::cerr << "[Error] 没有读取到任何采样点。" << std::endl;
        return -1;
    }

    // ===== 4) 确保用于保存路径文件的目录存在 =====
    std::string pathDir = "path";
    if (!fs::exists(pathDir))
    {
        if (!fs::create_directory(pathDir))
        {
            std::cerr << "[Error] 无法创建目录: " << pathDir << std::endl;
            return -1;
        }
    }

    // ===== 5) 打开一个日志文件 =====
    std::string reportFileName = pathDir + "/plan_report.txt";
    std::ofstream logFile(reportFileName);
    if (!logFile.is_open())
    {
        std::cerr << "[Error] 无法打开日志文件: " << reportFileName << std::endl;
        return -1;
    }
    logFile << "=== VFRRT Planner Report ===\n"
            << "[Info] 程序开始运行\n";

    // ===== 6) 设置规划器的参数 =====
    
    double exploration = 0.2;
    std::cout << "exploration: " << exploration << std::endl;
    double initial_lambda = 250.0;
    std::cout << "initial_lambda: " << initial_lambda << std::endl;
    unsigned int initial_lambda_samples = 1000000;
    std::cout << "initial_lambda_samples: " << initial_lambda_samples << std::endl;
    // ===== 7) 限制并发线程数 =====
    unsigned int maxConcurrentThreads = std::thread::hardware_concurrency();
    if (maxConcurrentThreads == 0) maxConcurrentThreads = 4;

    std::cout << "使用的最大并发线程数: " << maxConcurrentThreads << std::endl;
    logFile << "[Info] 使用的最大并发线程数: " << maxConcurrentThreads << "\n";

    // 用于管理所有的 future
    std::vector<std::future<PathPlanningResult>> futures;
    futures.reserve(samplePoints.size());

    // 用于保存所有的规划结果
    std::vector<PathPlanningResult> allResults;
    allResults.reserve(samplePoints.size());

    // 互斥锁，用于同步输出 (终端 & 日志)
    std::mutex coutMutex;

    // ===== 8) 提交多线程任务 =====
    int pathIndex = 0;
    for (const auto &point : samplePoints)
    {
        // 如果当前已提交的任务数 >= maxConcurrentThreads，则先收集已完成的
        while (futures.size() >= maxConcurrentThreads)
        {
            bool gotAny = false;
            for (auto it = futures.begin(); it != futures.end(); )
            {
                // 查看此 future 是否完成
                if (it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready)
                {
                    try
                    {
                        // 取得结果，加入 allResults
                        allResults.push_back(it->get());
                    }
                    catch (const std::exception &e)
                    {
                        std::lock_guard<std::mutex> lock(coutMutex);
                        std::cerr << "[Exception] 线程异常: " << e.what() << std::endl;
                        logFile << "[Error] 线程异常: " << e.what() << "\n";

                        // 如果出现异常，也存一个默认结果
                        PathPlanningResult dummy;
                        dummy.foundSolution = false;
                        dummy.pathIndex = -1;
                        dummy.isExactSolution = false;
                        allResults.push_back(dummy);
                    }
                    it = futures.erase(it);
                    gotAny = true;
                }
                else
                {
                    ++it;
                }
            }
            if (!gotAny)
            {
                // 若没有任何线程完成，就稍微睡一下，避免忙等待
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        }

        // 提交一个新任务
        futures.emplace_back(std::async(
            std::launch::async,
            [&vectorField, &point, &pathDir, pathIndex,
             exploration, initial_lambda, initial_lambda_samples,
             &u_combined, &v_combined, &coutMutex, &logFile]()
            {
                return planAndSavePath(
                    vectorField, point, pathDir, pathIndex,
                    exploration, initial_lambda, initial_lambda_samples,
                    u_combined, v_combined, coutMutex, logFile);
            }
        ));
        pathIndex++;
    }

    // ===== 9) 提交完全部任务后，收集剩余的结果 =====
    for (auto &fut : futures)
    {
        try
        {
            allResults.push_back(fut.get());
        }
        catch (const std::exception &e)
        {
            std::lock_guard<std::mutex> lock(coutMutex);
            std::cerr << "[Exception] 线程异常: " << e.what() << std::endl;
            logFile << "[Error] 线程异常: " << e.what() << "\n";

            PathPlanningResult dummy;
            dummy.foundSolution = false;
            dummy.pathIndex = -1;
            dummy.isExactSolution = false;
            allResults.push_back(dummy);
        }
    }
    futures.clear();

    // ===== 10) 从 allResults 中仅寻找 "精确解" 中上游成本最小的路径 =====
    double minExactCost = std::numeric_limits<double>::infinity();
    int minExactCostIndex = -1;

    for (const auto &res : allResults)
    {
        // 只比较精确解 (isExactSolution == true) 的成本
        if (res.foundSolution && res.isExactSolution && res.upstreamCost < minExactCost)
        {
            minExactCost = res.upstreamCost;
            minExactCostIndex = res.pathIndex;
        }
    }

    // 统计程序运行时间
    auto programEnd = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = programEnd - programStart;

    // ===== 输出最终结果 =====
    if (minExactCostIndex >= 0)
    {
        std::cout << "所有路径规划完成，最小上游成本的【精确解】路径索引: "
                  << minExactCostIndex
                  << "，其成本: " << minExactCost << std::endl;
        logFile << "[Result] 所有路径规划完成，最小上游成本的【精确解】路径索引: "
                << minExactCostIndex << "，成本: " << minExactCost << "\n";
    }
    else
    {
        std::cout << "所有路径规划完成，但没有找到任何【精确解】。" << std::endl;
        logFile << "[Result] 所有路径规划完成，但没有找到任何【精确解】。\n";
    }

    std::cout << "程序运行时间: " << elapsedSeconds.count() << " 秒" << std::endl;
    logFile << "[Info] 程序运行时间: " << elapsedSeconds.count() << " 秒\n";

    logFile.close();
    return 0;
}

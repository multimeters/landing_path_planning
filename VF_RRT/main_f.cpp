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

namespace ob = ompl::base;
namespace og = ompl::geometric;

using VectorField = std::function<Eigen::VectorXd(const ob::State *)>;

Eigen::MatrixXd loadMatrixFromFile(const std::string &filename, int rows, int cols)
{
    Eigen::MatrixXd matrix(rows, cols);
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return matrix;
    }

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            file >> matrix(i, j);
        }
    }

    file.close();
    return matrix;
}

void savePathToFile(const og::PathGeometric &path, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file for writing." << std::endl;
        return;
    }

    for (size_t i = 0; i < path.getStateCount(); ++i)
    {
        const auto *state = path.getState(i)->as<ob::RealVectorStateSpace::StateType>();
        file << state->values[0] << " " << state->values[1] << std::endl;
    }

    file.close();
}
void saveVectorFieldToFile(const VectorField &vectorField, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file for writing." << std::endl;
        return;
    }

    const int numPoints = 20; // Number of sample points along each axis
    const double step = 100.0 / (numPoints - 1);  // Step size for sampling points

    auto space = std::make_shared<ob::RealVectorStateSpace>(2);
    ob::ScopedState<ob::RealVectorStateSpace> state(space);

    for (double x = -50; x <= 50; x += step)
    {
        for (double y = -50; y <= 50; y += step)
        {
            state->values[0] = x;
            state->values[1] = y;

            Eigen::VectorXd field = vectorField(state.get());

            file << x << " " << y << " " << field[0] << " " << field[1] << std::endl;
        }
    }

    file.close();
}
int main()
{
    // Define the state space
    auto space = std::make_shared<ob::RealVectorStateSpace>(2);

    // Set the bounds for the R^2 space
    ob::RealVectorBounds bounds(2);
    bounds.setLow(-50);
    bounds.setHigh(50);
    space->setBounds(bounds);

    // Create a SimpleSetup object
    og::SimpleSetup ss(space);

    // Set the state validity checker (all states are valid in this simple example)
    ss.setStateValidityChecker([](const ob::State *state) { return true; });

    // Define the start and goal states
    ob::ScopedState<ob::RealVectorStateSpace> start(space);
    start->values[0] = -20;
    start->values[1] = -20;

    ob::ScopedState<ob::RealVectorStateSpace> goal(space);
    goal->values[0] = 20;
    goal->values[1] = 20;

    // Load vector field data from files
    Eigen::MatrixXd u_subset = loadMatrixFromFile("u_subset.txt", 4096, 2048);
    Eigen::MatrixXd v_subset = loadMatrixFromFile("v_subset.txt", 4096, 2048);
    Eigen::MatrixXd grad_x_norm = loadMatrixFromFile("grad_x_norm.txt", 4096, 2048);
    Eigen::MatrixXd grad_y_norm = loadMatrixFromFile("grad_y_norm.txt", 4096, 2048);

    // Define the vector field using the loaded data
    VectorField vectorField = [&u_subset, &v_subset, &grad_x_norm, &grad_y_norm](const ob::State *state) -> Eigen::VectorXd {
        const auto *realState = state->as<ob::RealVectorStateSpace::StateType>();
        double x = realState->values[0];
        double y = realState->values[1];

        // Map the coordinates to matrix indices (assuming -50 to 50 maps to matrix indices)
        int i = static_cast<int>(std::round((y + 50) / 100.0 * 4095)); // Row index
        int j = static_cast<int>(std::round((x + 50) / 100.0 * 2047)); // Column index

        // Ensure indices are within bounds
        i = std::min(std::max(i, 0), 4095);
        j = std::min(std::max(j, 0), 2047);

        Eigen::VectorXd field(2);
        field[0] = u_subset(i, j) + grad_x_norm(i, j);
        field[1] = v_subset(i, j) + grad_y_norm(i, j);
        return field;
    };
    // Save the vector field to a file
    saveVectorFieldToFile(vectorField, "vector_field.txt");
    // Set the parameters for VFRRT
    double exploration = 0.7; // Example value
    double initial_lambda = 17.5; // Example value
    unsigned int initial_lambda_samples = 100; // Example value

    // Create the VFRRT planner and set the vector field using its constructor
    auto planner = std::make_shared<og::VFRRT>(ss.getSpaceInformation(), vectorField, exploration, initial_lambda, initial_lambda_samples);

    // Set the planner for the SimpleSetup object
    ss.setPlanner(planner);

    // Set the start and goal states
    ss.setStartAndGoalStates(start, goal);

    // Attempt to solve the problem within a given time
    ob::PlannerStatus solved = ss.solve(10.0);

    if (solved)
    {
        std::cout << "Found solution:" << std::endl;
        ss.getSolutionPath().print(std::cout);

        // Save the solution path to a file
        savePathToFile(ss.getSolutionPath(), "path.txt");
    }
    else
    {
        std::cout << "No solution found" << std::endl;
    }

    return 0;
}

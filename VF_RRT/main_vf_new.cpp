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

    const int numPoints = 1024;  // Number of sample points along each axis
    const double step = 4096.0 / (numPoints - 1);  // Step size based on matrix dimensions

    auto space = std::make_shared<ob::RealVectorStateSpace>(2);
    ob::ScopedState<ob::RealVectorStateSpace> state(space);

    for (double x = 0; x <= 4096; x += step)
    {
        for (double y = 0; y <= 2048; y += step)
        {
            state->values[0] = x * 0.5;  // Convert matrix index to meters
            state->values[1] = y * 0.5;  // Convert matrix index to meters

            Eigen::VectorXd field = vectorField(state.get());

            file << x * 0.5 << " " << y * 0.5 << " " << field[0] << " " << field[1] << std::endl;
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
    bounds.setLow(-1024);  // Corresponding to -50 meters
    bounds.setHigh(1024);  // Corresponding to 50 meters
    space->setBounds(bounds);

    // Create a SimpleSetup object
    og::SimpleSetup ss(space);

    // Set the state validity checker (all states are valid in this simple example)
    ss.setStateValidityChecker([](const ob::State *state) { return true; });

    // Define the start and goal states
    ob::ScopedState<ob::RealVectorStateSpace> start(space);
    start->values[0] = 500;  // Corresponding to -20 meters in index space
    start->values[1] = 800;  // Corresponding to -20 meters in index space

    ob::ScopedState<ob::RealVectorStateSpace> goal(space);
    goal->values[0] = 1000;  // Corresponding to 20 meters in index space
    goal->values[1] = 400;  // Corresponding to 20 meters in index space

    // Load vector field data from files
    Eigen::MatrixXd u_combined = loadMatrixFromFile("/home/lhl/Amphi-RRT/VF_RRT/build/combined_u_vector_field_downsampled.txt", 4096, 2048);
    Eigen::MatrixXd v_combined = loadMatrixFromFile("/home/lhl/Amphi-RRT/VF_RRT/build/combined_v_vector_field_downsampled.txt", 4096, 2048);

    // Define the vector field using the loaded data
    VectorField vectorField = [&u_combined, &v_combined](const ob::State *state) -> Eigen::VectorXd {
        const auto *realState = state->as<ob::RealVectorStateSpace::StateType>();
        double x = realState->values[0] / 0.5;  // Convert meters to matrix index
        double y = realState->values[1] / 0.5;  // Convert meters to matrix index

        // Map the coordinates to matrix indices
        int i = static_cast<int>(std::round(y)); // Row index
        int j = static_cast<int>(std::round(x)); // Column index

        // Ensure indices are within bounds
        i = std::min(std::max(i, 0), 4095);
        j = std::min(std::max(j, 0), 2047);

        Eigen::VectorXd field(2);
        field[0] = u_combined(i, j);
        field[1] = v_combined(i, j);
        return field;
    };

    // Save the vector field to a file
    saveVectorFieldToFile(vectorField, "vector_field.txt");

    // Set the parameters for VFRRT
    double exploration = 1000;  // Example value
    double initial_lambda = 157.5;  // Example value
    unsigned int initial_lambda_samples = 10000;  // Example value

    // Create the VFRRT planner and set the vector field using its constructor
    auto planner = std::make_shared<og::VFRRT>(ss.getSpaceInformation(), vectorField, exploration, initial_lambda, initial_lambda_samples);
     
    // Set the planner for the SimpleSetup object
    ss.setPlanner(planner);

    // Set the start and goal states
    ss.setStartAndGoalStates(start, goal);

    // Attempt to solve the problem within a given time
    ob::PlannerStatus solved = ss.solve(1000.0);

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

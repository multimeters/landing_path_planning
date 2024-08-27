#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/VFRRT.h>
#include <ompl/util/RandomNumbers.h>
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <functional>

namespace ob = ompl::base;
namespace og = ompl::geometric;

using VectorField = std::function<Eigen::VectorXd(const ob::State *)>;

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

    // Define a simple vector field as an example
    VectorField vectorField = [](const ob::State *state) -> Eigen::VectorXd {
        const auto *realState = state->as<ob::RealVectorStateSpace::StateType>();
        Eigen::VectorXd field(2);
        field[0] = -realState->values[1];  // Example: rotate 90 degrees clockwise
        field[1] = realState->values[0];
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

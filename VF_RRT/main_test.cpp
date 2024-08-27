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
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>

namespace ob = ompl::base;

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

int main()
{
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

    return 0;
}

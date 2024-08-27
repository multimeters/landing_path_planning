% Read vector field from file
vectorFieldData = load('vector_field.txt');
x = vectorFieldData(:, 1);
y = vectorFieldData(:, 2);
u = vectorFieldData(:, 3);
v = vectorFieldData(:, 4);

% Read path from file
pathData = load('path.txt');
pathX = pathData(:, 1);
pathY = pathData(:, 2);

% Create a quiver plot for the vector field
figure;
quiver(x, y, u, v);
hold on;

% Plot the path
plot(pathX, pathY, 'r', 'LineWidth', 2);
scatter(pathX(1), pathY(1), 'go', 'filled'); % Start point
scatter(pathX(end), pathY(end), 'bo', 'filled'); % Goal point

% Add labels and title
xlabel('X');
ylabel('Y');
title('Vector Field and Path');
legend('Vector Field', 'Path', 'Start', 'Goal');

hold off;

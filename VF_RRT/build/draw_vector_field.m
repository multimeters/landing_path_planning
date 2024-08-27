% Load the vector field data
data = load('vector_field.txt');

% Extract the x, y coordinates and the corresponding vectors
x = data(:, 1);
y = data(:, 2);
vx = data(:, 3);
vy = data(:, 4);

% Plot the vector field
figure;
quiver(x, y, vx, vy);
axis equal;
xlabel('X');
ylabel('Y');
title('Vector Field');

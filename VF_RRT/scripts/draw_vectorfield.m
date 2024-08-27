% Load the data files
load('/home/lhl/Amphi-RRT/VF_RRT/result_data/total_sum_matrix.mat'); % Loads max_energy_matrix
load('/home/lhl/Amphi-RRT/VF_RRT/result_data/main_directions.mat'); % Loads main_directions
load('/home/lhl/Amphi-RRT/VF_RRT/result_data/dep_shoal_inlet_new_4096x2048.mat'); % Loads depth data

% Define the rows and columns of interest
rows_of_interest = 1000:3000;
cols_of_interest = 800:1800;

% Extract the relevant portion of the data
max_energy_matrix_subset = total_sum_matrix(rows_of_interest, cols_of_interest);
main_directions_subset = main_directions(rows_of_interest, cols_of_interest);
depth_subset = data(rows_of_interest, cols_of_interest); % Assuming depth is the variable name in the loaded file

% Compute the inverse (reciprocal) of the magnitudes
inverse_max_energy = 1 ./ max_energy_matrix_subset;

% Handle potential infinity values if any element in max_energy_matrix_subset is zero
inverse_max_energy(isinf(inverse_max_energy)) = 0;

% Normalize the inverse magnitudes
max_value = max(inverse_max_energy(:));
if (max_value > 0)
    normalized_inverse_energy = inverse_max_energy / max_value;
else
    normalized_inverse_energy = inverse_max_energy; % Handle the case when max_value is 0
end

% Convert directions to Cartesian components for the subset
u_subset = normalized_inverse_energy .* cos(main_directions_subset); % X-component (u)
v_subset = normalized_inverse_energy .* sin(main_directions_subset); % Y-component (v)

% Compute gradients of depth where the depth is less than 0.75
[grad_x, grad_y] = gradient(depth_subset);
depth_mask = depth_subset < 0.75;

% Use the gradients only where the depth is less than 0.75
grad_x(~depth_mask) = 0;
grad_y(~depth_mask) = 0;

% Normalize the gradients
magnitude = sqrt(grad_x.^2 + grad_y.^2);
grad_x_norm = grad_x ./ magnitude;
grad_y_norm = grad_y ./ magnitude;

% Handle cases where magnitude is zero to avoid NaNs in the normalization
grad_x_norm(magnitude == 0) = 0;
grad_y_norm(magnitude == 0) = 0;

% Create a grid for plotting the subset and gradients
[X_subset, Y_subset] = meshgrid(cols_of_interest, rows_of_interest);

% Define the file paths for saving
u_file = '/home/lhl/Amphi-RRT/VF_RRT/result_data/u_subset.txt';
v_file = '/home/lhl/Amphi-RRT/VF_RRT/result_data/v_subset.txt';
grad_x_file = '/home/lhl/Amphi-RRT/VF_RRT/result_data/grad_x_norm.txt';
grad_y_file = '/home/lhl/Amphi-RRT/VF_RRT/result_data/grad_y_norm.txt';

% Save each matrix to a text file using space as the delimiter
dlmwrite(u_file, u_subset, 'delimiter', ' ', 'precision', '%.6f');
dlmwrite(v_file, v_subset, 'delimiter', ' ', 'precision', '%.6f');
dlmwrite(grad_x_file, grad_x_norm, 'delimiter', ' ', 'precision', '%.6f');
dlmwrite(grad_y_file, grad_y_norm, 'delimiter', ' ', 'precision', '%.6f');


% Plot the combined vector fields on the same graph
figure;
hold on;
quiver(X_subset, Y_subset, u_subset, v_subset, 'b', 'AutoScaleFactor', 0.5); % Plotting the reciprocal magnitudes in blue
quiver(X_subset, Y_subset, grad_x_norm, grad_y_norm, 'r', 'AutoScaleFactor', 0.5); % Plotting the normalized gradients in red
axis equal;
title('Combined Vector Fields (Rows 1000 to 3000, Columns 800 to 1800)');
xlabel('X');
ylabel('Y');
legend('Reciprocal Magnitude Vector Field', 'Normalized Gradient Vector Field');
hold off;

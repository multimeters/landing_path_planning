% Load the data
data = load('dep_shoal_inlet.txt');

% Extract the first 256 rows
data_subset = data(1:256, :);

% Create a smooth gradient from -10 to 10 for each row
n_columns = size(data_subset, 2);
gradient = linspace(-10, 10, n_columns);  % Create a smooth gradient from -10 to 10

% Apply this gradient to each row, ensuring a smooth transition across columns
for i = 1:size(data_subset, 1)
    data_subset(i, :) = gradient;  % Assign the smooth gradient to each row
end

% Optionally apply Gaussian smoothing for further smoothness (if desired)
smoothed_data = imgaussfilt(data_subset, 1);  % You can adjust the filter size (1) for more or less smoothing

% Save the resulting data to a new text file
fileID = fopen('smoothed_depth_data.txt', 'w');  % Open a file for writing
for i = 1:size(smoothed_data, 1)
    fprintf(fileID, '%f\t', smoothed_data(i, :));  % Write each row with tab separation
    fprintf(fileID, '\n');  % New line after each row
end
fclose(fileID);  % Close the file

% Optionally, create a heatmap to visualize the result
figure;
imagesc(smoothed_data);  % Use imagesc to display the matrix as a heatmap
colorbar;                % Add a color bar to show the depth scale
title('Smoothed Depth Heatmap from Coast to Sea');
xlabel('Column Index');
ylabel('Row Index');
colormap('jet');         % You can choose a different colormap if needed

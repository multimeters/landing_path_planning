% Load the data from the file
data = load('dep_shoal_inlet.txt');

% Extract the first 256 rows
data_subset = data(1:256, :);

% Set all values from the 432nd column onward to -10.0
data_subset(:, 432:end) = -10.0;

% Swap the columns (reverse the order of the columns)
data_subset = data_subset(:, end:-1:1);

% Save the resulting matrix to a new text file without scientific notation
fileID = fopen('dep_shoal_inlet_subset_swapped.txt', 'w');
for i = 1:size(data_subset, 1)
    fprintf(fileID, '%f\t', data_subset(i, :));  % %f ensures no scientific notation
    fprintf(fileID, '\n');  % New line after each row
end
fclose(fileID);

% Create a heatmap of the extracted data
figure;
imagesc(data_subset);  % Use imagesc to display the matrix as a heatmap
colorbar;              % Add a color bar to show the depth scale
title('Depth Heatmap of the First 256 Rows with Swapped Columns');
xlabel('Column Index');
ylabel('Row Index');
colormap('jet');       % You can choose a different colormap if needed

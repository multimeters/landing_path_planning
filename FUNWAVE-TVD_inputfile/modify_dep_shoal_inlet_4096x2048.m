% 指定文件名
filename = 'dep_shoal_inlet_new.txt';

% 读取数据
data = -1*importdata(filename);

% 检查数据大小
[rows, cols] = size(data);
fprintf('数据大小为 %d 行 %d 列\n', rows, cols);

% 确保数据读取正确
if rows ~= 1024 || cols ~= 512
    error('数据尺寸不符合预期，请检查数据文件和读取过程。');
end

% 生成插值后新的行和列坐标
new_rows = 4096;
new_cols = 2048;
[x, y] = meshgrid(1:cols, 1:rows);
[xq, yq] = meshgrid(linspace(1, cols, new_cols), linspace(1, rows, new_rows));

% 对数据进行插值
new_data = -interp2(x, y, data, xq, yq, 'linear');
new_data = new_data.*1.5;
% 保存新的数据到txt文件，并设置格式为 0.000
new_filename = 'dep_shoal_inlet_new_4096x2048.txt';
fid = fopen(new_filename, 'w');
for i = 1:new_rows
    fprintf(fid, '%.3f\t', new_data(i, :));
    fprintf(fid, '\n');
end
fclose(fid);

% 绘制插值后的曲面图
figure; % 创建一个新的图形窗口
surf(new_data);
title('高程曲面图 (插值后)');
xlabel('列号');
ylabel('行号');
zlabel('高程值');
colorbar; % 显示颜色条以便查看高程值的分布
axis equal;  % 设置坐标轴比例相等

fprintf('插值后的数据已保存到文件 %s\n', new_filename);

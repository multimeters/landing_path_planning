% MATLAB 程序: visualize_path_and_vector_field_with_depth.m
% 该程序读取向量场、路径和深度数据文件，并进行可视化

% 清理环境
clear; clc; close all;

%% 1. 读取向量场数据
vectorFieldFile = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/vector_Field.txt';

% 检查文件是否存在
if ~isfile(vectorFieldFile)
    error('向量场文件 "%s" 不存在。请确保文件路径正确。', vectorFieldFile);
end

% 读取向量场数据
% 假设文件格式为: x y u v 每行一个采样点
vectorData = load(vectorFieldFile);

% 检查数据维度
if size(vectorData, 2) ~= 4
    error('向量场文件应包含四列数据: x y u v。');
end

x_vec = vectorData(:,1); % X 坐标（米）
y_vec = vectorData(:,2); % Y 坐标（米）
u = vectorData(:,3);     % 向量场在 X 方向的分量
v = vectorData(:,4);     % 向量场在 Y 方向的分量

% 根据分辨率确定唯一的 x 和 y 值
unique_x = unique(x_vec);
unique_y = unique(y_vec);
num_x = length(unique_x);
num_y = length(unique_y);

% 验证数据是否符合预期的网格尺寸
expected_num_x = 512; % X 轴采样点数
expected_num_y = 250; % Y 轴采样点数

if num_x ~= expected_num_x || num_y ~= expected_num_y
    warning('向量场的采样点数与预期不符。检查分辨率或数据是否正确。');
end

% 重新排列 u 和 v 为矩阵形式
% 根据 C++ 代码的保存顺序，外层循环是 Y，内层循环是 X
% 因此，在 MATLAB 中需要先按 Y，再按 X 排列
U = reshape(u, [num_x, num_y])';
V = reshape(v, [num_x, num_y])';
X = reshape(x_vec, [num_x, num_y])';
Y = reshape(y_vec, [num_x, num_y])';

%% 2. 读取路径数据
pathFile = '\\wsl.localhost\Ubuntu\/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/path/path_2.txt';

% 检查文件是否存在
if ~isfile(pathFile)
    error('路径文件 "%s" 不存在。请确保文件路径正确。', pathFile);
end

% 读取路径数据
% 假设文件格式为: x y 每行一个路径点
pathData = load(pathFile);

% 检查数据维度
if size(pathData, 2) ~= 2
    error('路径文件应包含两列数据: x y。');
end

path_x = pathData(:,1); % 路径点的 X 坐标（米）
path_y = pathData(:,2); % 路径点的 Y 坐标（米）

% 确定起点和终点
start_x = path_x(1);
start_y = path_y(1);
end_x = path_x(end);
end_y = path_y(end);

%% 3. 读取深度数据
depthFile = '\\wsl.localhost\ubuntu\home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/depth_a15_512x512.txt';

% 检查文件是否存在
if ~isfile(depthFile)
    error('深度文件 "%s" 不存在。请确保文件路径正确。', depthFile);
end

% 读取深度数据
% 假设 depth_a15.txt 是一个 250x512 的二维数组，按行排列 Y 方向，按列排列 X 方向
% 每个元素代表对应网格点的深度值（米）
Depth = load(depthFile);

% 检查深度数据维度
if ~ismatrix(Depth) || size(Depth,1) ~= 256 || size(Depth,2) ~= 512
    error('深度文件应为一个 250x512 的二维数组。');
end

% 定义分辨率
dx = 1; % X 方向分辨率（米）
dy = 2; % Y 方向分辨率（米）

% 生成对应的 X 和 Y 网格
[n, m] = size(Depth);
xx = [0:m-1] * dx; % X 轴坐标
yy = [0:n-1] * dy; % Y 轴坐标

% 创建 X 和 Y 网格用于绘图
[X_depth, Y_depth] = meshgrid(xx, yy);

%% 4. 可视化

% 设置图形窗口属性
wid = 6; % 宽度（英寸）
len = 6; % 高度（英寸）
figure;
set(gcf, 'Units', 'inches', ...
         'PaperUnits', 'inches', ...
         'PaperSize', [wid len], ...
         'Position', [1 1 wid len], ...
         'PaperPosition', [0 0 wid len]);
hold on;
axis equal;
grid on;

% 绘制深度热力图
% 使用 pcolor 绘制深度热力图，并应用 shading flat
% 为了与向量场和路径对齐，使用 -Depth 使深度为负值
p = pcolor(X_depth, Y_depth, -Depth);
shading flat;
colormap('jet');
colorbar;
caxis([min(-Depth(:)) max(-Depth(:))]); % 根据实际深度范围调整

% 绘制深度等高线
% 绘制多个深度等高线，例如从 -10 到 0，每隔1米一条
contour_levels = -10:1:0; % 深度范围，根据实际数据调整
[~, h_contour] = contour(X_depth, Y_depth, -Depth, contour_levels, 'w');
set(h_contour, 'LineWidth', 0.5); % 设置等高线线宽

% 特别标注深度小于 0.75 米的等高线
% 深度 < 0.75 米 对应 -Depth > -0.75
contour_specific = -0.75;
[~, h_specific] = contour(X_depth, Y_depth, -Depth, [contour_specific contour_specific], 'w--', 'LineWidth', 2, 'DisplayName', '深度 < 0.75 米');

% 绘制向量场
% 为了避免向量过于密集，可以选择间隔绘制
% 例如，每隔20个点绘制一次
skip_x = 5; % 每隔20个X点绘制一次
skip_y = 5; % 每隔20个Y点绘制一次

quiver(X(1:skip_y:end, 1:skip_x:end), ...
       Y(1:skip_y:end, 1:skip_x:end), ...
       U(1:skip_y:end, 1:skip_x:end), ...
       V(1:skip_y:end, 1:skip_x:end), ...
       'AutoScale', 'on', 'AutoScaleFactor', 0.4, 'Color', [0.7 0.7 0.7], ...
       'MaxHeadSize', 0.2, 'LineWidth', 0.2, 'DisplayName', '向量场');

% 绘制路径
plot(path_x, path_y, 'b-', 'LineWidth', 2, 'DisplayName', '路径');

% 标记起点
plot(start_x, start_y, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', '起点');

% 标记终点
plot(end_x, end_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', '终点');

% 添加图例
legend('深度热力图', '深度等高线', '深度 < 0.75 米', '向量场', '路径', '起点', '终点', 'Location', 'best');

% 添加标题和标签
title('路径规划、向量场与深度热力图可视化');
xlabel('X 位置 (米)');
ylabel('Y 位置 (米)');

hold off;

% 可选：保存图像
% saveas(gcf, 'path_vector_depth_visualization.png');

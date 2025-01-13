% MATLAB 程序: visualize_all_paths_with_vector_field_and_depth.m
% 该程序读取向量场、所有路径和深度数据文件，并进行可视化
% 修改：绘制所有路径，其中最优路径 path37 用绿色并加粗绘制，其终点为最佳登陆点
% 新增：绘制采样的海岸线点，并确保最佳登陆点绘制在采样点上方

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

%% 2. 读取所有路径数据
pathFolder = '\\wsl.localhost\Ubuntu\home\lhl\share\rip\scripts\landing_path_planning\VFRRT_Planner\path\';
pathPattern = fullfile(pathFolder, 'path_*.txt');
pathFiles = dir(pathPattern);

% 检查是否找到路径文件
if isempty(pathFiles)
    error('路径文件夹 "%s" 中未找到符合 "path_*.txt" 格式的文件。', pathFolder);
end

% 识别最优路径文件（假设为 path_37.txt）
optimalPathFileName = 'path_37.txt';
isOptimalPathAvailable = false;
for i = 1:length(pathFiles)
    if strcmp(pathFiles(i).name, optimalPathFileName)
        isOptimalPathAvailable = true;
        break;
    end
end

%% 3. 读取深度数据
depthFile = '\\wsl.localhost\ubuntu\home\lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/depth_a15_512x512.txt';

% 检查文件是否存在
if ~isfile(depthFile)
    error('深度文件 "%s" 不存在。请确保文件路径正确。', depthFile);
end

% 读取深度数据
% 假设 depth_a15.txt 是一个 256x512 的二维数组，按行排列 Y 方向，按列排列 X 方向
% 每个元素代表对应网格点的深度值（米）
Depth = load(depthFile);

% 检查深度数据维度
expected_depth_rows = 256; % 根据注释应为256x512
expected_depth_cols = 512;
if ~ismatrix(Depth) || size(Depth,1) ~= expected_depth_rows || size(Depth,2) ~= expected_depth_cols
    error('深度文件应为一个 %dx%d 的二维数组。', expected_depth_rows, expected_depth_cols);
end

% 定义分辨率
dx = 1; % X 方向分辨率（米）
dy = 2; % Y 方向分辨率（米）

% 生成对应的 X 和 Y 网格
[n, m] = size(Depth);
xx = (0:m-1) * dx; % X 轴坐标
yy = (0:n-1) * dy; % Y 轴坐标

% 创建 X 和 Y 网格用于绘图
[X_depth, Y_depth] = meshgrid(xx, yy);

%% 4. 读取并准备绘制采样的海岸线数据
% 新增部分：读取采样的海岸线数据
coastlineFile = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/sampled_coastline.txt';

% 检查文件是否存在
if ~isfile(coastlineFile)
    warning('海岸线采样文件 "%s" 不存在。将跳过绘制海岸线。', coastlineFile);
    sampledCoastline = [];
else
    % 读取海岸线数据
    % 假设文件格式为: x y 每行一个采样点
    sampledCoastline = load(coastlineFile);
    
    % 检查数据维度
    if size(sampledCoastline, 2) ~= 2
        warning('海岸线采样文件 "%s" 应包含两列数据: x y。将跳过绘制海岸线。', coastlineFile);
        sampledCoastline = [];
    end
end

%% 5. 可视化

% 设置图形窗口属性
wid = 14; % 宽度（英寸）
len = 10; % 高度（英寸）
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
% 更改等高线颜色为黑色，以便在图例中可见
contour_levels = -10:1:0; % 深度范围，根据实际数据调整
[~, h_contour] = contour(X_depth, Y_depth, -Depth, contour_levels, 'k');
set(h_contour, 'LineWidth', 0.5); % 设置等高线线宽

% 特别标注深度小于 0.75 米的等高线
% 深度 < 0.75 米 对应 -Depth > -0.75
contour_specific = -0.75;
[~, h_specific] = contour(X_depth, Y_depth, -Depth, [contour_specific contour_specific], 'k--', 'LineWidth', 2, 'DisplayName', '深度 < 0.75 米');

% 绘制向量场
% 为了避免向量过于密集，可以选择间隔绘制
% 例如，每隔5个点绘制一次
skip_x = 5; % 每隔5个X点绘制一次
skip_y = 5; % 每隔5个Y点绘制一次

quiver_handle = quiver(X(1:skip_y:end, 1:skip_x:end), ...
           Y(1:skip_y:end, 1:skip_x:end), ...
           U(1:skip_y:end, 1:skip_x:end), ...
           V(1:skip_y:end, 1:skip_x:end), ...
           'AutoScale', 'on', 'AutoScaleFactor', 0.4, 'Color', [0.7 0.7 0.7], ...
           'MaxHeadSize', 0.2, 'LineWidth', 0.2, 'DisplayName', '向量场');

% 绘制采样的海岸线点（绘制在路径之前，以便路径和登陆点覆盖在其上）
if ~isempty(sampledCoastline)
    plot(sampledCoastline(:,1), sampledCoastline(:,2), 'ro', 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'DisplayName', '采样海岸线点');
end

% 绘制所有路径
% 创建一个颜色映射，用于为不同路径分配不同颜色
numPaths = length(pathFiles);
colorMap = lines(numPaths); % 使用 MATLAB 的 lines 颜色图

legend_entries = {}; % 初始化图例条目
legend_handles = []; % 初始化图例句柄

for i = 1:numPaths
    currentPathFile = pathFiles(i).name;
    currentPathFullPath = fullfile(pathFolder, currentPathFile);
    
    % 读取路径数据
    pathData = load(currentPathFullPath);
    
    % 检查数据维度
    if size(pathData, 2) ~= 2
        warning('路径文件 "%s" 应包含两列数据: x y。跳过此文件。', currentPathFile);
        continue;
    end
    
    pathX = pathData(:,1);
    pathY = pathData(:,2);
    
    if strcmp(currentPathFile, optimalPathFileName)
        % 绘制最优路径（绿色加粗）
        h = plot(pathX, pathY, 'g-', 'LineWidth', 2.5, 'DisplayName', '最优路径 (path37)');
        
        % 标记最优路径的终点为最佳登陆点（绘制在采样点之后，确保在其上方）
        best_landing_x = pathX(end);
        best_landing_y = pathY(end);
        plot(best_landing_x, best_landing_y, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', '最佳登陆点');
        
        % 标记起点
        start_x = pathX(1);
        start_y = pathY(1);
        plot(start_x, start_y, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', '起点');
        
        % 添加到图例句柄和条目
        legend_handles = [legend_handles, h];
        legend_entries = [legend_entries, '最优路径 (path37)'];
        legend_entries = [legend_entries, '最佳登陆点', '起点'];
    else
        % 绘制其他路径
        h = plot(pathX, pathY, 'Color', colorMap(i,:), 'LineWidth', 1, 'DisplayName', sprintf('路径 %d', i));
        
        % 添加到图例句柄和条目（只为第一条普通路径添加图例，以避免图例过长）
        if i == 1
            legend_handles = [legend_handles, h];
            legend_entries = [legend_entries, '其他路径'];
        end
    end
end

% 绘制最优路径已经在循环中完成，无需重复

% 绘制采样的海岸线点已经在之前完成

% 添加图例
% 手动添加图例条目，避免重复和过多条目
legend_labels = {'深度热力图', '等深线', '深度 < 0.75 米', '向量场', '采样海岸线点', '其他路径', '最优路径 (path37)', '最佳登陆点', '起点'};

% 获取图例句柄
h_pcolor = p; % 深度热力图
h_contour_main = h_contour; % 深度等高线
h_contour_specific = h_specific; % 特别标注深度等高线
h_quiver = quiver_handle; % 向量场
if ~isempty(sampledCoastline)
    h_coastline = findobj(gca, 'DisplayName', '采样海岸线点');
else
    h_coastline = [];
end
h_otherPaths = findobj(gca, 'DisplayName', '其他路径');
h_optimalPath = findobj(gca, 'DisplayName', '最优路径 (path37)');
h_bestLanding = findobj(gca, 'DisplayName', '最佳登陆点');
h_start = findobj(gca, 'DisplayName', '起点');

% 创建图例句柄数组
legend_handles = [h_pcolor, h_contour_main, h_contour_specific, h_quiver];
if ~isempty(sampledCoastline)
    legend_handles = [legend_handles, h_coastline];
end
legend_handles = [legend_handles, h_otherPaths, h_optimalPath, h_bestLanding, h_start];

% 创建图例标签数组
legend_labels = {'深度热力图', '等深线', '深度 < 0.75 米', '向量场'};
if ~isempty(sampledCoastline)
    legend_labels = [legend_labels, '采样海岸线点'];
end
legend_labels = [legend_labels, '其他路径', '最优路径 (path37)', '最佳登陆点', '起点'];

legend(legend_handles, legend_labels, 'Location', 'best');

% 添加标题和标签
title('路径规划、向量场与深度热力图可视化');
xlabel('X 位置 (米)');
ylabel('Y 位置 (米)');

hold off;

% 可选：保存图像
% saveas(gcf, 'all_paths_vector_depth_visualization.png'); % 如果需要保存图像，请取消注释

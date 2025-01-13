% MATLAB程序：绘制路径及每个点的切向量，并计算与x轴的夹角（以弧度表示）
% 同时将x和y坐标转换为索引（x对应列索引，y对应行索引），
% 并与main_directions中的方向计算差异，保存相关数据到指定文件夹

% 清空工作区和图形窗口
clear;
clc;
close all;

% 指定输入文件名
input_filename = 'sampled_path.txt';

% 指定输出文件夹名称
output_folder = 'output_data'; % 你可以根据需要更改文件夹名称

% 指定输出文件名
angles_filename = 'tangent_angles.txt';
direction_differences_filename = 'direction_differences.txt'; % 方向差异文件名
figure_filename = 'path_with_tangents.png'; % 可选：图形保存的文件名

% 创建输出文件夹（如果不存在）
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('已创建文件夹 %s。\n', output_folder);
else
    fprintf('文件夹 %s 已存在。\n', output_folder);
end

% 构建完整的输入文件路径
input_filepath = fullfile(pwd, input_filename);

% 检查输入文件是否存在
if ~isfile(input_filepath)
    error('文件 %s 不存在。请确保文件在当前工作目录中。', input_filename);
end

% 读取数据
% 假设文件中的数据是数值型，并且每行包含两个数值（x和y坐标）
try
    data = load(input_filepath); % 使用load函数读取数据
catch
    error('无法读取文件 %s。请确保文件格式正确。', input_filename);
end

% 检查数据维度
if size(data, 2) < 2
    error('文件 %s 的数据格式不正确。至少需要两列数据（x和y坐标）。', input_filename);
end

% 提取x和y坐标
x = data(:,1);
y = data(:,2);

% 检查点的数量
num_points = length(x);
if num_points < 2
    error('路径至少需要两个点来计算切向量。当前点数：%d', num_points);
end

% 初始化切向量
tangent_x = zeros(num_points, 1);
tangent_y = zeros(num_points, 1);

% 计算切向量
for i = 1:num_points
    if i == 1
        % 第一个点，使用前向差分
        dx = x(2) - x(1);
        dy = y(2) - y(1);
    elseif i == num_points
        % 最后一个点，使用后向差分
        dx = x(end) - x(end-1);
        dy = y(end) - y(end-1);
    else
        % 中间点，使用中心差分
        dx = x(i+1) - x(i-1);
        dy = y(i+1) - y(i-1);
    end
    tangent_x(i) = dx;
    tangent_y(i) = dy;
end

% 设置箭头的基准长度
vector_length = 20; % 根据需要调整

% 归一化切向量并缩放
tangent_norm = sqrt(tangent_x.^2 + tangent_y.^2);
% 避免除以零
tangent_norm(tangent_norm == 0) = 1;
tangent_unit_x = (tangent_x ./ tangent_norm) * vector_length;
tangent_unit_y = (tangent_y ./ tangent_norm) * vector_length;

% 计算每个切向量与x轴的夹角（以弧度表示）
% 使用atan2(y, x)函数，结果范围为 (-pi, pi)
angles_rad = atan2(tangent_y, tangent_x);

% 将角度存储到输出文件夹中的文件
angles_filepath = fullfile(output_folder, angles_filename);
try
    % 打开文件以写入（如果文件存在，将被覆盖）
    fileID = fopen(angles_filepath, 'w');
    if fileID == -1
        error('无法打开文件 %s 进行写入。', angles_filepath);
    end

    % 写入每个角度，格式为每行一个角度，固定小数点格式
    for i = 1:num_points
        fprintf(fileID, '%.6f\n', angles_rad(i));
    end

    % 关闭文件
    fclose(fileID);
    fprintf('切向量与x轴的夹角已成功写入文件 %s。\n', angles_filepath);
catch ME
    error('写入文件时出错: %s', ME.message);
end

% 将位置和角度保存到一个方向差异文件中（使用固定小数点格式）
direction_differences_filepath = fullfile(output_folder, direction_differences_filename);
try
    % 将x和y转换为索引
    index_x = floor(x / 1);  % x除以1取整数，对应列索引
    index_y = floor(y / 2);  % y除以2取整数，对应行索引

    % 将索引转换为1-based，以适应MATLAB的索引方式
    index_x_mapped = index_x + 1; % 列索引
    index_y_mapped = index_y + 1; % 行索引

    % 获取main_directions变量路径
    main_directions_path = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/main_directions.mat';

    % 检查main_directions.mat文件是否存在
    if ~isfile(main_directions_path)
        error('文件 %s 不存在。请确保文件路径正确。', main_directions_path);
    end

    % 加载main_directions
    try
        main_data = load(main_directions_path);
    catch
        error('无法加载文件 %s。请确保文件格式正确。', main_directions_path);
    end

    % 假设main_directions变量存在于main_data中
    if ~isfield(main_data, 'main_directions')
        error('文件 %s 中不存在变量 main_directions。', main_directions_path);
    end

    main_directions = main_data.main_directions;

    % 获取main_directions的尺寸
    [max_index_y, max_index_x] = size(main_directions); % 行数对应y方向，列数对应x方向

    % 创建一个逻辑掩码，标记哪些点的索引在main_directions范围内
    valid_mask = (index_x_mapped >= 1) & (index_x_mapped <= max_index_x) & ...
                 (index_y_mapped >= 1) & (index_y_mapped <= max_index_y);

    num_valid = sum(valid_mask);
    num_invalid = num_points - num_valid;

    fprintf('共有 %d 个点在main_directions范围内，%d 个点超出范围将被忽略。\n', num_valid, num_invalid);

    if num_valid == 0
        error('没有任何点的索引在main_directions的范围内。请检查数据和main_directions的尺寸。');
    end

    % 过滤有效点的数据
    index_x_valid = index_x(valid_mask);
    index_y_valid = index_y(valid_mask);
    angles_rad_valid = angles_rad(valid_mask);
    index_x_mapped_valid = index_x_mapped(valid_mask);
    index_y_mapped_valid = index_y_mapped(valid_mask);

    % 初始化main_direction和delta_direc
    main_direction = zeros(num_valid, 1);
    delta_direc = zeros(num_valid, 1);

    % 计算main_direction和delta_direc
    for i = 1:num_valid
        % 获取当前点的索引
        ix = index_x_mapped_valid(i); % 列索引
        iy = index_y_mapped_valid(i); % 行索引

        % 获取main_direction
        main_direction_i = main_directions(iy, ix);

        % 存储main_direction
        main_direction(i) = main_direction_i;

        % 计算delta_direc = angles_rad - main_direction
        delta = angles_rad_valid(i) - main_direction_i;

        % 将delta_direc规范化到[-pi, pi]
        delta = mod(delta + pi, 2*pi) - pi;

        % 存储delta_direc
        delta_direc(i) = abs(delta);
    end

    % 打开方向差异文件以写入（如果文件存在，将被覆盖）
    fileID_delta = fopen(direction_differences_filepath, 'w');
    if fileID_delta == -1
        error('无法打开文件 %s 进行写入。', direction_differences_filepath);
    end

    % 写入每个有效点的 index_x, index_y, angles_rad, main_direction, delta_direc
    for i = 1:num_valid
        fprintf(fileID_delta, '%d %d %.6f %.6f %.6f\n', ...
            index_x_valid(i), index_y_valid(i), angles_rad_valid(i), main_direction(i), delta_direc(i));
    end

    % 关闭方向差异文件
    fclose(fileID_delta);
    fprintf('方向差异已成功保存到文件 %s。\n', direction_differences_filepath);
catch ME
    error('保存方向差异文件时出错: %s', ME.message);
end

% 设置 quiver 的缩放参数
scale_quiver = 0.5; % 根据需要调整

% 绘制路径
figure;
plot(x, y, '-b', 'LineWidth', 2); % 用蓝色线条绘制路径
hold on;

% 绘制所有切向量
quiver(x, y, tangent_unit_x, tangent_unit_y, scale_quiver, 'g', ...
       'LineWidth', 1.5, 'MaxHeadSize', 1.5);

% 绘制所有点
plot(x, y, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k'); % 用黑色圆点标记所有点

% 高亮第一个点
first_x = x(1);
first_y = y(1);
plot(first_x, first_y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % 用红色圆点标记第一个点

% 添加图形标题和标签
title('路径及每个点的切向量');
xlabel('X 坐标');
ylabel('Y 坐标');

% 添加图例
legend('路径', '切向量', '所有点', '第一个点', 'Location', 'best');

% 添加网格
grid on;

% 设置轴比例相等，以防止箭头失真
axis equal;

% 保持图形
hold off;

% 保存图形到输出文件夹（可选）
figure_filepath = fullfile(output_folder, figure_filename);
saveas(gcf, figure_filepath);
fprintf('图形已保存到文件 %s。\n', figure_filepath);

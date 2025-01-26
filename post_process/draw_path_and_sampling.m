% process_all_paths.m
%
% 功能:
% 处理指定文件夹中的所有 path_*.txt 文件，包括路径采样、切向量计算、
% 方向差异分析及结果的保存和可视化。

%% ---------------------------- 初始化部分 ----------------------------

% 清空工作区和图形窗口
clear;
clc;
close all;

% --------------------- 输入和输出路径设置 ---------------------

% 定义包含所有 path 文件的文件夹路径
input_folder = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/path_e_05_lam_0011_sam_1000000/';

% 定义采样点输出的基础文件夹路径
sampled_output_base_folder = '\\wsl.localhost\ubuntu\home\lhl\share\rip\scripts\landing_path_planning\post_process\sampled_paths\';

% 定义数据处理和结果输出的基础文件夹路径
processed_output_base_folder = '\\wsl.localhost\ubuntu\home\lhl\share\rip\scripts\landing_path_planning\post_process\processed_data\';

% 定义 main_directions.mat 文件的路径
main_directions_path = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/main_directions.mat';

% 创建输出文件夹（如果不存在）
if ~exist(sampled_output_base_folder, 'dir')
    mkdir(sampled_output_base_folder);
    fprintf('已创建采样点输出文件夹 %s。\n', sampled_output_base_folder);
else
    fprintf('采样点输出文件夹 %s 已存在。\n', sampled_output_base_folder);
end

if ~exist(processed_output_base_folder, 'dir')
    mkdir(processed_output_base_folder);
    fprintf('已创建处理数据输出文件夹 %s。\n', processed_output_base_folder);
else
    fprintf('处理数据输出文件夹 %s 已存在。\n', processed_output_base_folder);
end

% --------------------- 加载 main_directions ---------------------

% 检查 main_directions.mat 文件是否存在
if ~isfile(main_directions_path)
    error('文件 %s 不存在。请确保文件路径正确。', main_directions_path);
end

% 加载 main_directions
try
    main_data = load(main_directions_path);
catch
    error('无法加载文件 %s。请确保文件格式正确。', main_directions_path);
end

% 检查 main_directions 变量是否存在
if ~isfield(main_data, 'main_directions')
    error('文件 %s 中不存在变量 main_directions。', main_directions_path);
end

main_directions = main_data.main_directions;

% 获取 main_directions 的尺寸
[max_main_y, max_main_x] = size(main_directions); % 行数对应 y 方向，列数对应 x 方向

%% ---------------------------- 处理所有 path 文件 ----------------------------

% 获取所有 path_*.txt 文件的列表
path_files = dir(fullfile(input_folder, 'path_*.txt'));

% 检查是否找到任何 path 文件
if isempty(path_files)
    error('在文件夹 %s 中未找到任何 path_*.txt 文件。', input_folder);
end

fprintf('在文件夹 %s 中找到 %d 个 path 文件。\n', input_folder, length(path_files));

% 循环处理每个 path 文件
for pf = 1:length(path_files)
    % 获取当前 path 文件的信息
    current_path_file = path_files(pf).name;
    current_path_full = fullfile(input_folder, current_path_file);
    
    fprintf('\n===== 处理文件 %s (%d/%d) =====\n', current_path_file, pf, length(path_files)));
    
    %% ---------------------------- 第一步：路径采样 ----------------------------
    
    % 读取数据
    % 假设文件中每行包含两个或三个数值，分别对应 x、y（和 z）坐标
    try
        data_original = load(current_path_full);
    catch
        warning('无法读取文件 %s。跳过该文件。', current_path_file);
        continue; % 跳过当前循环，继续下一个文件
    end
    
    % 检查数据格式
    [numRows_orig, numCols_orig] = size(data_original);
    if numCols_orig < 2
        warning('数据格式错误：文件 %s 应至少包含两列，分别表示 X 和 Y 坐标。跳过该文件。', current_path_file);
        continue;
    end
    
    % 提取 X 和 Y（以及 Z）坐标
    x_orig = data_original(:, 1);
    y_orig = data_original(:, 2);
    if numCols_orig >= 3
        z_orig = data_original(:, 3);
    else
        z_orig = [];
    end
    
    % 绘制原始路径
    figure('Name', sprintf('原始路径 - %s', current_path_file));
    plot(x_orig, y_orig, '-b', 'LineWidth', 2, 'DisplayName', '原始路径');
    hold on;
    
    xlabel('X 坐标');
    ylabel('Y 坐标');
    title(sprintf('原始路径规划线 - %s', current_path_file));
    grid on;
    axis equal; % 保持坐标轴比例一致
    
    % 如果存在 Z 坐标，绘制三维图形
    if ~isempty(z_orig)
        figure('Name', sprintf('原始路径三维图 - %s', current_path_file));
        plot3(x_orig, y_orig, z_orig, '-b', 'LineWidth', 2, 'DisplayName', '原始路径');
        xlabel('X 坐标');
        ylabel('Y 坐标');
        zlabel('Z 坐标');
        title(sprintf('三维路径规划线 - %s', current_path_file));
        grid on;
        axis equal;
    end
    
    % 沿路径每隔1米进行采样
    fprintf('开始沿路径每隔1米进行采样...\n');
    
    % 计算每段的距离
    if isempty(z_orig)
        % 2D 路径
        deltas_orig = diff([x_orig, y_orig], 1, 1);
        segment_lengths_orig = sqrt(deltas_orig(:,1).^2 + deltas_orig(:,2).^2);
    else
        % 3D 路径
        deltas_orig = diff([x_orig, y_orig, z_orig], 1, 1);
        segment_lengths_orig = sqrt(deltas_orig(:,1).^2 + deltas_orig(:,2).^2 + deltas_orig(:,3).^2);
    end
    
    % 计算累计距离
    cumulative_distance_orig = [0; cumsum(segment_lengths_orig)];
    
    % 总路径长度
    total_length_orig = cumulative_distance_orig(end);
    
    % 定义采样间隔
    sample_interval = 1; % 米
    
    % 定义采样点的位置
    sample_distances = 0:sample_interval:total_length_orig;
    
    % 如果最后一个采样点不在总长度上，则添加总长度
    if sample_distances(end) < total_length_orig
        sample_distances = [sample_distances, total_length_orig];
    end
    
    % 进行插值以找到采样点的坐标
    if isempty(z_orig)
        % 2D 路径
        sampled_x = interp1(cumulative_distance_orig, x_orig, sample_distances);
        sampled_y = interp1(cumulative_distance_orig, y_orig, sample_distances);
        sampled_points = [sampled_x', sampled_y'];
    else
        % 3D 路径
        sampled_x = interp1(cumulative_distance_orig, x_orig, sample_distances);
        sampled_y = interp1(cumulative_distance_orig, y_orig, sample_distances);
        sampled_z = interp1(cumulative_distance_orig, z_orig, sample_distances);
        sampled_points = [sampled_x', sampled_y', sampled_z'];
    end
    
    % 绘制采样点
    if isempty(z_orig)
        plot(sampled_x, sampled_y, 'r*', 'MarkerSize', 8, 'DisplayName', '采样点');
    else
        plot(sampled_x, sampled_y, 'r*', 'MarkerSize', 8, 'DisplayName', '采样点');
        figure(1); % 原始路径 2D 图
        hold on;
        plot3(sampled_x, sampled_y, sampled_z, 'r*', 'MarkerSize', 8, 'DisplayName', '采样点');
    end
    legend('原始路径', '采样点');
    hold off;
    
    % 保存采样点到文本文件
    % 根据当前 path 文件名生成唯一的采样输出文件名
    [~, path_name, ~] = fileparts(current_path_file); % 获取不带扩展名的文件名
    sampled_path_file = fullfile(sampled_output_base_folder, sprintf('%s_sampled.txt', path_name));
    
    try
        % 创建保存目录（如果不存在）
        [savePath_sampled, ~, ~] = fileparts(sampled_path_file);
        if ~exist(savePath_sampled, 'dir')
            mkdir(savePath_sampled);
        end
    
        % 保存格式：每行一个点，列分别为 X Y (Z)
        if isempty(z_orig)
            save(sampled_path_file, 'sampled_points', '-ascii');
        else
            save(sampled_path_file, 'sampled_points', '-ascii');
        end
    
        fprintf('采样完成。采样点已保存到: %s\n', sampled_path_file);
    catch ME
        warning('保存采样点到文件 %s 时出错: %s。', sampled_path_file, ME.message);
    end
    
    %% ---------------------------- 第二步：路径处理 ----------------------------
    
    % 读取采样后的路径文件
    input_filepath = sampled_path_file; % 使用采样后的路径文件
    
    % 检查输入文件是否存在
    if ~isfile(input_filepath)
        warning('文件 %s 不存在。跳过该文件。', input_filepath);
        continue;
    end
    
    % 读取数据
    % 假设文件中的数据是数值型，并且每行包含两个数值（x和y坐标）或三个（x,y,z）
    try
        data_sampled = load(input_filepath); % 使用 load 函数读取数据
    catch
        warning('无法读取文件 %s。跳过该文件。', input_filepath);
        continue;
    end
    
    % 检查数据维度
    [num_points, num_cols_sampled] = size(data_sampled);
    if num_cols_sampled < 2
        warning('文件 %s 的数据格式不正确。至少需要两列数据（x和y坐标）。跳过该文件。', input_filepath);
        continue;
    end
    
    % 提取 x 和 y 坐标
    x = data_sampled(:,1);
    y = data_sampled(:,2);
    if num_cols_sampled >=3
        z = data_sampled(:,3);
    else
        z = [];
    end
    
    % 检查点的数量
    if num_points < 2
        warning('路径至少需要两个点来计算切向量。文件 %s 当前点数：%d。跳过该文件。', input_filepath, num_points);
        continue;
    end
    
    % 初始化切向量
    tangent_x = zeros(num_points, 1);
    tangent_y = zeros(num_points, 1);
    
    % 计算切向量
    fprintf('计算切向量...\n');
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
    % 使用 atan2(y, x) 函数，结果范围为 (-pi, pi)
    angles_rad = atan2(tangent_y, tangent_x);
    
    % 将角度存储到输出文件夹中的文件
    angles_filename = sprintf('%s_tangent_angles.txt', path_name);
    angles_filepath = fullfile(processed_output_base_folder, angles_filename);
    
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
        warning('写入文件 %s 时出错: %s', angles_filepath, ME.message);
    end
    
    %% --------------------- 方向差异计算 ---------------------
    
    fprintf('计算方向差异...\n');
    
    % 将 x 和 y 转换为索引
    % 根据实际情况调整缩放因子（例如，x 对应列索引，y 对应行索引）
    scale_x = 1; % x 除以 scale_x 取整数，对应列索引
    scale_y = 2; % y 除以 scale_y 取整数，对应行索引
    
    index_x = floor(x / scale_x);  % 列索引
    index_y = floor(y / scale_y);  % 行索引
    
    % 将索引转换为 1-based，以适应 MATLAB 的索引方式
    index_x_mapped = index_x + 1; % 列索引
    index_y_mapped = index_y + 1; % 行索引
    
    % 创建一个逻辑掩码，标记哪些点的索引在 main_directions 范围内
    valid_mask = (index_x_mapped >= 1) & (index_x_mapped <= max_main_x) & ...
                 (index_y_mapped >= 1) & (index_y_mapped <= max_main_y);
    
    num_valid = sum(valid_mask);
    num_invalid = num_points - num_valid;
    
    fprintf('共有 %d 个点在 main_directions 范围内，%d 个点超出范围将被忽略。\n', num_valid, num_invalid);
    
    if num_valid == 0
        warning('没有任何点的索引在 main_directions 的范围内。跳过方向差异计算。');
        continue;
    end
    
    % 过滤有效点的数据
    index_x_valid = index_x(valid_mask);
    index_y_valid = index_y(valid_mask);
    angles_rad_valid = angles_rad(valid_mask);
    index_x_mapped_valid = index_x_mapped(valid_mask);
    index_y_mapped_valid = index_y_mapped(valid_mask);
    
    % 初始化 main_direction 和 delta_direc
    main_direction = zeros(num_valid, 1);
    delta_direc = zeros(num_valid, 1);
    
    % 计算 main_direction 和 delta_direc
    for i = 1:num_valid
        % 获取当前点的索引
        ix = index_x_mapped_valid(i); % 列索引
        iy = index_y_mapped_valid(i); % 行索引
    
        % 获取 main_direction
        main_direction_i = main_directions(iy, ix);
    
        % 存储 main_direction
        main_direction(i) = main_direction_i;
    
        % 计算 delta_direc = angles_rad - main_direction
        delta = angles_rad_valid(i) - main_direction_i;
    
        % 将 delta_direc 规范化到 [-pi, pi]
        delta = mod(delta + pi, 2*pi) - pi;
    
        % 存储 delta_direc
        delta_direc(i) = abs(delta);
    end
    
    % 将方向差异保存到文件
    direction_differences_filename = sprintf('%s_direction_differences.txt', path_name);
    direction_differences_filepath = fullfile(processed_output_base_folder, direction_differences_filename);
    
    try
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
        warning('保存方向差异文件 %s 时出错: %s', direction_differences_filepath, ME.message);
    end
    
    %% ---------------------------- 第三步：绘制路径及切向量 ----------------------------
    
    fprintf('绘制路径及切向量...\n');
    
    % 设置 quiver 的缩放参数
    scale_quiver = 0.5; % 根据需要调整
    
    % 绘制路径
    figure('Name', sprintf('路径及切向量 - %s', current_path_file));
    plot(x, y, '-b', 'LineWidth', 2, 'DisplayName', '路径');
    hold on;
    
    % 绘制所有切向量
    quiver(x, y, tangent_unit_x, tangent_unit_y, scale_quiver, 'g', ...
           'LineWidth', 1.5, 'MaxHeadSize', 1.5, 'DisplayName', '切向量');
    
    % 绘制所有点
    plot(x, y, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'DisplayName', '所有点'); % 用黑色圆点标记所有点
    
    % 高亮第一个点
    first_x = x(1);
    first_y = y(1);
    plot(first_x, first_y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '第一个点'); % 用红色圆点标记第一个点
    
    % 添加图形标题和标签
    title(sprintf('路径及每个点的切向量 - %s', current_path_file));
    xlabel('X 坐标');
    ylabel('Y 坐标');
    
    % 添加图例
    legend('Location', 'best');
    
    % 添加网格
    grid on;
    
    % 设置轴比例相等，以防止箭头失真
    axis equal;
    
    % 保持图形
    hold off;
    
    % 保存图形到输出文件夹
    figure_filename = sprintf('%s_path_with_tangents.png', path_name);
    figure_filepath = fullfile(processed_output_base_folder, figure_filename);
    
    try
        saveas(gcf, figure_filepath);  
        fprintf('图形已保存到文件 %s。\n', figure_filepath);
    catch ME
        warning('无法保存图形到文件 %s: %s', figure_filepath, ME.message);
    end
    
    % 关闭所有图形窗口以节省内存
    close all;
    
    fprintf('===== 完成文件 %s 的处理 =====\n', current_path_file);
end

fprintf('\n===== 所有文件处理完毕 =====\n');

function combined_plot_and_process_path(original_path_file, sampled_path_file, main_directions_path, varargin)
% combined_plot_and_process_path
% 
% 功能:
% 1. 读取指定路径文件，绘制路径线，并沿路径每隔指定间隔采样并保存采样点。
% 2. 读取采样点，绘制路径及每个点的切向量，计算与x轴的夹角（以弧度表示）。
%    同时将x和y坐标转换为索引，并与main_directions中的方向计算差异，
%    保存相关数据到指定文件夹。
%
% 输入参数:
%   original_path_file     - 原始路径文件的完整路径
%   sampled_path_file      - 采样点输出文件的完整路径
%   main_directions_path   - main_directions.mat 文件的完整路径
%
% 可选参数（通过名称-值对指定）:
%   'OutputFolder'                - 输出文件夹的名称或路径 (默认: 'output_data')
%   'SampleInterval'              - 采样间隔，单位为米 (默认: 1)
%   'FigureFilename'              - 保存绘图的文件名 (默认: 'path_with_tangents.png')
%   'AnglesFilename'              - 保存切向量夹角的文件名 (默认: 'tangent_angles.txt')
%   'DirectionDifferencesFilename' - 保存方向差异的文件名 (默认: 'direction_differences.txt')
%   'ScaleX'                      - x坐标转换为索引的缩放因子 (默认: 1)
%   'ScaleY'                      - y坐标转换为索引的缩放因子 (默认: 2)
%
% 示例:
%   combined_plot_and_process_path(original_path, sampled_path, main_directions, ...
%       'OutputFolder', 'results', 'SampleInterval', 0.5, ...
%       'AnglesFilename', 'angles_path0.txt', ...
%       'DirectionDifferencesFilename', 'differences_path0.txt');

    %% ---------------------------- 输入参数解析 ----------------------------

    p = inputParser;
    addRequired(p, 'original_path_file', @ischar);
    addRequired(p, 'sampled_path_file', @ischar);
    addRequired(p, 'main_directions_path', @ischar);
    addParameter(p, 'OutputFolder', 'output_data', @ischar);
    addParameter(p, 'SampleInterval', 1, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'FigureFilename', 'path_with_tangents.png', @ischar);
    addParameter(p, 'AnglesFilename', 'tangent_angles.txt', @ischar);
    addParameter(p, 'DirectionDifferencesFilename', 'direction_differences.txt', @ischar);
    addParameter(p, 'ScaleX', 1, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'ScaleY', 2, @(x) isnumeric(x) && x > 0);

    parse(p, original_path_file, sampled_path_file, main_directions_path, varargin{:});

    output_folder = p.Results.OutputFolder;
    sample_interval = p.Results.SampleInterval;
    figure_filename = p.Results.FigureFilename;
    angles_filename = p.Results.AnglesFilename;
    direction_differences_filename = p.Results.DirectionDifferencesFilename;
    scale_x = p.Results.ScaleX;
    scale_y = p.Results.ScaleY;

    %% ---------------------------- 初始化部分 ----------------------------

    % 清空不必要的变量，保留输入参数
    clearvars -except original_path_file sampled_path_file main_directions_path output_folder ...
               sample_interval figure_filename angles_filename direction_differences_filename ...
               scale_x scale_y;
    clc;
    close all;

    %% --------------------- 创建输出文件夹 ---------------------

    % 创建输出文件夹（如果不存在）
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
        fprintf('已创建文件夹 %s。\n', output_folder);
    else
        fprintf('文件夹 %s 已存在。\n', output_folder);
    end

    %% ---------------------------- 第一步：路径采样 ----------------------------

    fprintf('\n===== 第一步：路径采样 =====\n');

    % 检查原始路径文件是否存在
    if ~isfile(original_path_file)
        error('指定的文件不存在: %s', original_path_file);
    end

    % 读取数据
    % 假设文件中每行包含两个或三个数值，分别对应 x、y（和 z）坐标
    try
        data_original = load(original_path_file);
    catch
        error('无法读取文件 %s。请确保文件格式正确。', original_path_file);
    end

    % 检查数据格式
    [numRows_orig, numCols_orig] = size(data_original);
    if numCols_orig < 2
        error('数据格式错误：文件应至少包含两列，分别表示 X 和 Y 坐标。');
    end

    % 提取 X 和 Y（以及 Z）坐标
    x_orig = data_original(:, 1);
    y_orig = data_original(:, 2);
    if numCols_orig >= 3
        z_orig = data_original(:, 3);
    else
        z_orig = [];
    end

    % **新增：移除连续重复的点**
    if isempty(z_orig)
        % 2D 路径
        duplicate_mask = [false; (x_orig(2:end) == x_orig(1:end-1)) & (y_orig(2:end) == y_orig(1:end-1))];
    else
        % 3D 路径
        duplicate_mask = [false; (x_orig(2:end) == x_orig(1:end-1)) & ...
                              (y_orig(2:end) == y_orig(1:end-1)) & ...
                              (z_orig(2:end) == z_orig(1:end-1))];
    end

    if any(duplicate_mask)
        fprintf('检测到 %d 个连续重复的点，已移除这些点。\n', sum(duplicate_mask));
        x_orig(duplicate_mask) = [];
        y_orig(duplicate_mask) = [];
        if ~isempty(z_orig)
            z_orig(duplicate_mask) = [];
        end
    else
        fprintf('未检测到连续重复的点。\n');
    end

    % 绘制原始路径（图形不可见）
    fig1 = figure('Name', '原始路径', 'Visible', 'off');
    plot(x_orig, y_orig, '-b', 'LineWidth', 2, 'DisplayName', '原始路径');
    hold on;

    xlabel('X 坐标');
    ylabel('Y 坐标');
    title('原始路径规划线');
    grid on;
    axis equal; % 保持坐标轴比例一致

    % 如果存在 Z 坐标，绘制三维图形（图形不可见）
    if ~isempty(z_orig)
        fig2 = figure('Name', '原始路径三维图', 'Visible', 'off');
        plot3(x_orig, y_orig, z_orig, '-b', 'LineWidth', 2, 'DisplayName', '原始路径');
        xlabel('X 坐标');
        ylabel('Y 坐标');
        zlabel('Z 坐标');
        title('三维路径规划线');
        grid on;
        axis equal;
    end

    % 沿路径每隔 sample_interval 米进行采样
    fprintf('开始沿路径每隔 %.2f 米进行采样...\n', sample_interval);

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

    % 定义采样点的位置
    sample_distances = 0:sample_interval:total_length_orig;

    % 如果最后一个采样点不在总长度上，则添加总长度
    if sample_distances(end) < total_length_orig
        sample_distances = [sample_distances, total_length_orig];
    end

    % **新增：确保 cumulative_distance_orig 是严格递增的**
    % 通过移除重复的 cumulative_distance_orig 值
    [cumulative_distance_orig_unique, unique_idx] = unique(cumulative_distance_orig, 'stable');
    if length(cumulative_distance_orig_unique) ~= length(cumulative_distance_orig)
        warning('cumulative_distance_orig 中存在重复值，已移除重复项。');
        cumulative_distance_orig = cumulative_distance_orig_unique;
        if isempty(z_orig)
            x_orig = x_orig(unique_idx);
            y_orig = y_orig(unique_idx);
        else
            x_orig = x_orig(unique_idx);
            y_orig = y_orig(unique_idx);
            z_orig = z_orig(unique_idx);
        end
    end

    % 进行插值以找到采样点的坐标
    if isempty(z_orig)
        % 2D 路径
        sampled_x = interp1(cumulative_distance_orig, x_orig, sample_distances, 'linear');
        sampled_y = interp1(cumulative_distance_orig, y_orig, sample_distances, 'linear');
        sampled_points = [sampled_x', sampled_y'];
    else
        % 3D 路径
        sampled_x = interp1(cumulative_distance_orig, x_orig, sample_distances, 'linear');
        sampled_y = interp1(cumulative_distance_orig, y_orig, sample_distances, 'linear');
        sampled_z = interp1(cumulative_distance_orig, z_orig, sample_distances, 'linear');
        sampled_points = [sampled_x', sampled_y', sampled_z'];
    end

    % 绘制采样点（在原始路径图形中绘制）
    if isempty(z_orig)
        plot(sampled_x, sampled_y, 'r*', 'MarkerSize', 8, 'DisplayName', '采样点');
    else
        plot(sampled_x, sampled_y, 'r*', 'MarkerSize', 8, 'DisplayName', '采样点');
        figure(fig1);
        hold on;
        plot3(sampled_x, sampled_y, sampled_z, 'r*', 'MarkerSize', 8, 'DisplayName', '采样点');
    end
    legend('Original Path', 'Sampled Points');
    hold off;

    % 保存采样点到文本文件
    try
        % 创建保存目录（如果不存在）
        [savePath_sampled, ~, ~] = fileparts(sampled_path_file);
        if ~exist(savePath_sampled, 'dir')
            mkdir(savePath_sampled);
        end

        % 保存格式：每行一个点，列分别为 X Y (Z)
        save(sampled_path_file, 'sampled_points', '-ascii');

        fprintf('采样完成。采样点已保存到: %s\n', sampled_path_file);
    catch ME
        error('保存采样点时出错: %s', ME.message);
    end

    %% ---------------------------- 第二步：路径处理 ----------------------------

    fprintf('\n===== 第二步：路径处理 =====\n');

    % 读取采样后的路径文件
    input_filepath = sampled_path_file; % 使用采样后的路径文件

    % 检查输入文件是否存在
    if ~isfile(input_filepath)
        error('文件 %s 不存在。请确保文件路径正确。', input_filepath);
    end

    % 读取数据
    % 假设文件中的数据是数值型，并且每行包含两个数值（x和y坐标）或三个（x,y,z）
    try
        data_sampled = load(input_filepath); % 使用load函数读取数据
    catch
        error('无法读取文件 %s。请确保文件格式正确。', input_filepath);
    end

    % 检查数据维度
    [num_points, num_cols_sampled] = size(data_sampled);
    if num_cols_sampled < 2
        error('文件 %s 的数据格式不正确。至少需要两列数据（x和y坐标）。', input_filepath);
    end

    % 提取x和y坐标
    x = data_sampled(:,1);
    y = data_sampled(:,2);
    if num_cols_sampled >=3
        z = data_sampled(:,3);
    else
        z = [];
    end

    % 检查点的数量
    if num_points < 2
        error('路径至少需要两个点来计算切向量。当前点数：%d', num_points);
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

    % --------------------- 方向差异计算 ---------------------

    fprintf('计算方向差异...\n');

    % 将x和y转换为索引
    % 你可以根据实际情况调整缩放因子（例如，x对应列索引，y对应行索引）
    index_x = floor(x / scale_x);  % 列索引
    index_y = floor(y / scale_y);  % 行索引

    % 将索引转换为1-based，以适应MATLAB的索引方式
    index_x_mapped = index_x + 1; % 列索引
    index_y_mapped = index_y + 1; % 行索引

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

    % 将方向差异保存到文件
    direction_differences_filepath = fullfile(output_folder, direction_differences_filename);
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
        error('保存方向差异文件时出错: %s', ME.message);
    end

    %% ---------------------------- 第三步：绘制路径及切向量 ----------------------------

    fprintf('\n===== 第三步：绘制路径及切向量 =====\n');

    % 设置 quiver 的缩放参数
    scale_quiver = 0.5; % 根据需要调整

    % 绘制路径及切向量（图形不可见）
    fig3 = figure('Name', '路径及切向量', 'Visible', 'off');
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
    title('路径及每个点的切向量');
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

    % 保存图形到输出文件夹（可选）
    figure_filepath = fullfile(output_folder, figure_filename);
    try
        saveas(fig3, figure_filepath);  
        fprintf('图形已保存到文件 %s。\n', figure_filepath);
    catch ME
        warning('无法保存图形到文件 %s: %s', figure_filepath, ME.message);
    end

    % 关闭所有图形窗口以释放资源
    close(fig1);
    if ~isempty(z_orig)
        close(fig2);
    end
    close(fig3);

    fprintf('\n===== 脚本执行完毕 =====\n');
end

function visualize_all_paths_with_vector_field_and_depth(outputPath)
    % 该函数读取向量场、路径数据、深度数据、海岸线数据，并生成一张包含所有路径、向量场和深度信息的可视化图像
    % 参数：
    %   outputPath: 保存可视化图像的路径，例如 'C:/output/all_paths_visualization.png'

    % 清理环境
    %clear; clc; close all;

    %% 1. 读取向量场数据
    vectorFieldFile = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/vector_Field.txt';
    % 检查文件是否存在
    if ~isfile(vectorFieldFile)
        error('向量场文件 "%s" 不存在。请确保文件路径正确。', vectorFieldFile);
    end
    vectorData = load(vectorFieldFile);
    if size(vectorData, 2) ~= 4
        error('向量场文件应包含四列数据: x y u v。');
    end

    % 分离数据
    x_vec = vectorData(:,1);
    y_vec = vectorData(:,2);
    u = vectorData(:,3);
    v = vectorData(:,4);

    % 根据分辨率确定唯一的 x 和 y 值
    unique_x = unique(x_vec);
    unique_y = unique(y_vec);
    num_x = length(unique_x);
    num_y = length(unique_y);
    expected_num_x = 512;
    expected_num_y = 250;
    if num_x ~= expected_num_x || num_y ~= expected_num_y
        warning('向量场的采样点数与预期不符。检查分辨率或数据是否正确。');
    end

    % 重新排列 u 和 v 为矩阵形式
    U = reshape(u, [num_x, num_y])';
    V = reshape(v, [num_x, num_y])';
    X = reshape(x_vec, [num_x, num_y])';
    Y = reshape(y_vec, [num_x, num_y])';

    %% 2. 读取路径数据
    pathFolder = '\\wsl.localhost\Ubuntu\home\lhl\share\rip\scripts\landing_path_planning\VFRRT_Planner\path\';
    pathPattern = fullfile(pathFolder, 'path_*.txt');
    pathFiles = dir(pathPattern);
    if isempty(pathFiles)
        error('路径文件夹 "%s" 中未找到符合 "path_*.txt" 格式的文件。', pathFolder);
    end

    %% 3. 读取深度数据
    depthFile = '\\wsl.localhost\ubuntu\home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/depth_a15_512x512.txt';
    if ~isfile(depthFile)
        error('深度文件 "%s" 不存在。请确保文件路径正确。', depthFile);
    end
    Depth = load(depthFile);
    expected_depth_rows = 256;
    expected_depth_cols = 512;
    if size(Depth,1) ~= expected_depth_rows || size(Depth,2) ~= expected_depth_cols
        error('深度文件应为一个 %dx%d 的二维数组。', expected_depth_rows, expected_depth_cols);
    end
    dx = 1;
    dy = 2;
    [n, m] = size(Depth);
    xx = (0:m-1) * dx;
    yy = (0:n-1) * dy;
    [X_depth, Y_depth] = meshgrid(xx, yy);

    %% 4. 读取海岸线数据
    coastlineFile = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/data/sampled_coastline.txt';
    if isfile(coastlineFile)
        sampledCoastline = load(coastlineFile);
        if size(sampledCoastline, 2) ~= 2
            warning('海岸线采样文件 "%s" 应包含两列数据: x y。', coastlineFile);
            sampledCoastline = [];
        end
    else
        sampledCoastline = [];
    end

    %% 5. 可视化
    % 创建 figure 对象，并设置为不可见
    figure_handle = figure('Visible', 'off');
    hold on;
    axis equal;
    grid on;
    p = pcolor(X_depth, Y_depth, -Depth);
    shading flat;
    colormap('jet');
    colorbar;
    caxis([min(-Depth(:)) max(-Depth(:))]);
    contour_levels = -10:1:0;
    [~, h_contour] = contour(X_depth, Y_depth, -Depth, contour_levels, 'k');
    set(h_contour, 'LineWidth', 0.5);

    % 绘制向量场
    skip_x = 5;
    skip_y = 5;
    quiver_handle = quiver(X(1:skip_y:end, 1:skip_x:end), Y(1:skip_y:end, 1:skip_x:end), ...
           U(1:skip_y:end, 1:skip_x:end), V(1:skip_y:end, 1:skip_x:end), ...
           'AutoScale', 'on', 'AutoScaleFactor', 0.4, 'Color', [0.7 0.7 0.7], 'MaxHeadSize', 0.2);

    % 绘制海岸线点
    if ~isempty(sampledCoastline)
        plot(sampledCoastline(:,1), sampledCoastline(:,2), 'ro', 'MarkerSize', 1, 'MarkerFaceColor', 'r');
    end

    % 绘制路径
    colorMap = lines(length(pathFiles));
    for i = 1:length(pathFiles)
        currentPathFile = pathFiles(i).name;
        currentPathFullPath = fullfile(pathFolder, currentPathFile);
        pathData = load(currentPathFullPath);
        if size(pathData, 2) ~= 2
            warning('路径文件 "%s" 应包含两列数据: x y。');
            continue;
        end
        pathX = pathData(:,1);
        pathY = pathData(:,2);
        plot(pathX, pathY, 'Color', colorMap(i,:), 'LineWidth', 1);
    end

    % 添加图例
    legend({'深度热力图', '等深线', '向量场', '采样海岸线点', '路径'}, 'Location', 'best');
    title('路径规划与深度热力图');
    xlabel('X 位置 (米)');
    ylabel('Y 位置 (米)');

    % 保存图像
    saveas(figure_handle, outputPath);
    disp(['图像已保存到: ', outputPath]);

    % 关闭图形窗口
    close(figure_handle);
end

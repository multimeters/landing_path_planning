% 指定方向差异文件路径
directionFilePath = 'output_data/direction_differences_path_0.txt';

% 指定深度文件路径
depthFilePath = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/source/depth_a15_512x512.txt'; % 请将此路径替换为实际深度文件的路径

% 检查方向差异文件是否存在
if exist(directionFilePath, 'file')
    % 加载方向差异数据
    directionData = load(directionFilePath);
    
    % 获取尺寸
    directionSize = size(directionData);
    
    % 输出尺寸
    fprintf('文件 "%s" 的尺寸为：\n', directionFilePath);
    disp(directionSize);
    
    % 检查数据是否有至少5列
    if directionSize(2) >= 5
        % 提取第一行的第1、2和5个数据
        firstRow = directionData(1, :);
        element1 = firstRow(1); % 列索引
        element2 = firstRow(2); % 行索引
        element5 = firstRow(5); % μ
        
        % 显示提取的数据
        fprintf('第一行的第1个数据 (列索引): %d\n', element1);
        fprintf('第一行的第2个数据 (行索引): %d\n', element2);
        fprintf('第一行的第5个数据 (μ): %f\n', element5);
    else
        error('数据的列数不足5列。无法提取第5个数据。');
    end
else
    error('文件 "%s" 未找到。请检查路径是否正确。', directionFilePath);
end

% 检查深度文件是否存在
if exist(depthFilePath, 'file')
    % 加载深度数据
    depthData = load(depthFilePath);
    
    % 获取深度数据的尺寸
    depthSize = size(depthData);
    
    % 输出深度数据的尺寸
    fprintf('深度文件 "%s" 的尺寸为：\n', depthFilePath);
    disp(depthSize);
    
    % 验证 element1 和 element2 是否在深度数据的范围内
    if element2 >= 1 && element2 <= depthSize(1) && ...
       element1 >= 1 && element1 <= depthSize(2)
        % 提取对应的深度 h
        h = depthData(element2, element1);
        
        fprintf('提取的深度 h (element2=%d, element1=%d): %f 米\n', element2, element1, h);
    else
        error('element1 或 element2 超出了深度数据的范围。');
    end
else
    error('深度文件 "%s" 未找到。请检查路径是否正确。', depthFilePath);
end

% 参数定义
U_kmh = 10; % 速度 U，单位为 km/h
U = U_kmh * (1000/3600); % 将速度转换为 m/s

g = 9.81; % 重力加速度，单位为 m/s²

% 检查 h 是否为正数
if h <= 0
    error('深度 h 必须为正数。当前 h = %f', h);
end

% 假设 μ (element5) 是以度为单位，如果是弧度请跳过此步骤
mu_degrees = element5;
mu = deg2rad(mu_degrees); % 将 μ 转换为弧度

% 计算表达式 1 - (U * cos(mu)) / sqrt(g * h)
term = (U * cos(mu)) / sqrt(g * h);
result = 1 - term;

% 显示结果
fprintf('计算结果 1 - (U * cos(mu)) / sqrt(g * h) = %f\n', result);

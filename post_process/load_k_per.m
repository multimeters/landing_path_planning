% 指定方向差异文件路径
directionFilePath = 'output_data/direction_differences_path_0.txt';

% 指定深度文件路径
depthFilePath = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/source/depth_a15_512x512.txt'; % 请将此路径替换为实际深度文件的路径

% 指定 omega.mat 文件路径
omegaMatPath = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/source/omega.mat';

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

% ------------------- 新增部分 -------------------

% 检查 omega.mat 文件是否存在
if exist(omegaMatPath, 'file')
    % 加载 omega.mat 文件
    omegaData = load(omegaMatPath);
    
    % 假设 omega 的变量名为 'omega_1024'，如果不同，请相应更改
    if isfield(omegaData, 'omega_1024')
        omega = omegaData.omega_1024;
        fprintf('成功加载 omega 数组，尺寸为：\n');
        disp(size(omega));
    else
        error('omega.mat 中未找到变量 "omega_1024"。请检查变量名。');
    end
else
    error('omega.mat 文件 "%s" 未找到。请检查路径是否正确。', omegaMatPath);
end

% 检查 result 是否为数值且非空
if ~isnumeric(result) || isempty(result)
    error('变量 "result" 无效。请确保之前的计算成功。');
end

% 计算修正频率
% omega 是以弧度表示的频率，result 是标量
omega_corrected = omega * result;

% 检查 omega_corrected 是否有效（非零，避免除以零）
if any(omega_corrected == 0)
    warning('存在 omega_corrected 为零的元素，转换为周期时将导致无穷大。');
end

% 将修正频率转换为周期 T = 2 * pi / omega_corrected
T = (2 * pi) ./ omega_corrected;

% 显示部分修正频率和对应周期
numElementsToDisplay = min(10, length(T)); % 显示前10个元素
fprintf('前 %d 个修正频率和对应周期：\n', numElementsToDisplay);
for i = 1:numElementsToDisplay
    fprintf('omega_corrected(%d) = %f rad/s, T(%d) = %f s\n', i, omega_corrected(i), i, T(i));
end

% 处理 omega_corrected:
% 1. 每个元素除以0.1
% 2. 四舍五入取整
% 3. 如果结果为 Inf 或 -Inf，设为0

% 首先，创建一个复制的数组以进行处理
omega_processed = omega_corrected / 0.1;
mu_processed = mu / 0.5;
% 四舍五入取整
omega_processed = round(omega_processed);
mu_processed = round(mu_processed);
% 找出无限大的元素（Inf 或 -Inf）
inf_indices = isinf(omega_corrected);

% 将无限大的元素设为0
omega_processed(inf_indices) = 0;

% 显示部分处理后的结果
numElementsToDisplay = min(10, length(omega_processed)); % 显示前10个元素
fprintf('前 %d 个处理后的 omega_corrected 元素：\n', numElementsToDisplay);
for i = 1:numElementsToDisplay
    fprintf('omega_processed(%d) = %d\n', i, omega_processed(i));
end

% 如果需要，保存处理后的结果到文件
% save('omega_processed.mat', 'omega_processed');

% ------------------- 新增部分 -------------------

% 加载 heave.mat 文件
heaveMatPath = 'heave.mat'; % 请确保 heave.mat 的路径正确

if exist(heaveMatPath, 'file')
    % 加载 heave.mat 文件
    heaveData = load(heaveMatPath);
    
    % 假设 heave.mat 中的变量名为 'heave'，如果不同，请相应更改
    if isfield(heaveData, 'heave')
        heave = heaveData.heave;
        fprintf('成功加载 heave 数据，尺寸为：\n');
        disp(size(heave));
    else
        error('heave.mat 中未找到变量 "heave"。请检查变量名。');
    end
else
    error('heave.mat 文件 "%s" 未找到。请检查路径是否正确。', heaveMatPath);
end

% 检查 heave 数据是否至少有两列
if size(heave, 2) < 2
    error('heave 数据必须至少包含两列。');
end

% 初始化 heave_result 数组，与 omega_processed 具有相同长度
heave_result = zeros(length(omega_processed), 1);

% 遍历每个元素，查找匹配并赋值
for i = 1:length(omega_processed)
    % 当前的 omega 和 mu 值
    current_omega = omega_processed(i);
    current_mu = mu_processed;
    
    % 查找 heave 数据中第一列等于 current_omega 且第二列等于 current_mu 的行
    matches = (heave(:,1) == current_omega) & (heave(:,2) == current_mu);
    
    if any(matches)
        % 假设需要提取第三列的值，如果 heave 数据有更多列，可以根据需要调整
        heave_value = heave(matches, 4);
        
        % 如果有多个匹配，取第一个匹配的值
        heave_result(i) = heave_value(1);
    else
        % 如果找不到匹配，赋值为 0
        heave_result(i) = 0;
    end
end

% 显示部分 heave_result
numElementsToDisplay = min(10, length(heave_result)); % 显示前10个元素
fprintf('前 %d 个 heave_result 元素：\n', numElementsToDisplay);
for i = 1:numElementsToDisplay
    fprintf('heave_result(%d) = %d\n', i, heave_result(i));
end

% 如果需要，将 heave_result 保存到文件
% save('heave_result.mat', 'heave_result');

% ------------------- 结束新增部分 -------------------
% ------------------- 处理 roll.mat -------------------

% 加载 roll.mat 文件
rollMatPath = 'roll.mat'; % 请确保 roll.mat 的路径正确

if exist(rollMatPath, 'file')
    % 加载 roll.mat 文件
    rollDataStruct = load(rollMatPath);
    
    % 假设 roll.mat 中的变量名为 'roll'，如果不同，请相应更改
    if isfield(rollDataStruct, 'roll')
        roll = rollDataStruct.roll;
        fprintf('成功加载 roll 数据，尺寸为：\n');
        disp(size(roll));
    else
        error('roll.mat 中未找到变量 "roll"。请检查变量名。');
    end
else
    error('roll.mat 文件 "%s" 未找到。请检查路径是否正确。', rollMatPath);
end

% 检查 roll 数据是否至少有四列
if size(roll, 2) < 4
    error('roll 数据必须至少包含四列。');
end

% 初始化 roll_result 数组，与 omega_processed 具有相同长度
roll_result = zeros(length(omega_processed), 1);

% 遍历每个元素，查找匹配并赋值
for i = 1:length(omega_processed)
    % 当前的 omega 和 mu 值
    current_omega = omega_processed(i);
    current_mu = mu_processed;
    
    % 查找 roll 数据中第一列等于 current_omega 且第二列等于 current_mu 的行
    matches = (roll(:,1) == current_omega) & (roll(:,2) == current_mu);
    
    if any(matches)
        % 假设需要提取第四列的值，如果 roll 数据有更多列，可以根据需要调整
        roll_value = roll(matches, 4);
        
        % 如果有多个匹配，取第一个匹配的值
        roll_result(i) = roll_value(1);
    else
        % 如果找不到匹配，赋值为 0
        roll_result(i) = 0;
    end
end

% 显示部分 roll_result
numElementsToDisplay = min(10, length(roll_result)); % 显示前10个元素
fprintf('前 %d 个 roll_result 元素：\n', numElementsToDisplay);
for i = 1:numElementsToDisplay
    fprintf('roll_result(%d) = %d\n', i, roll_result(i));
end

% 如果需要，将 roll_result 保存到文件
% save('roll_result.mat', 'roll_result');

% ------------------- 处理 pitch.mat -------------------

% 加载 pitch.mat 文件
pitchMatPath = 'pitch.mat'; % 请确保 pitch.mat 的路径正确

if exist(pitchMatPath, 'file')
    % 加载 pitch.mat 文件
    pitchDataStruct = load(pitchMatPath);
    
    % 假设 pitch.mat 中的变量名为 'pitch'，如果不同，请相应更改
    if isfield(pitchDataStruct, 'pitch')
        pitch = pitchDataStruct.pitch;
        fprintf('成功加载 pitch 数据，尺寸为：\n');
        disp(size(pitch));
    else
        error('pitch.mat 中未找到变量 "pitch"。请检查变量名。');
    end
else
    error('pitch.mat 文件 "%s" 未找到。请检查路径是否正确。', pitchMatPath);
end

% 检查 pitch 数据是否至少有四列
if size(pitch, 2) < 4
    error('pitch 数据必须至少包含四列。');
end

% 初始化 pitch_result 数组，与 omega_processed 具有相同长度
pitch_result = zeros(length(omega_processed), 1);

% 遍历每个元素，查找匹配并赋值
for i = 1:length(omega_processed)
    % 当前的 omega 和 mu 值
    current_omega = omega_processed(i);
    current_mu = mu_processed;
    
    % 查找 pitch 数据中第一列等于 current_omega 且第二列等于 current_mu 的行
    matches = (pitch(:,1) == current_omega) & (pitch(:,2) == current_mu);
    
    if any(matches)
        % 假设需要提取第四列的值，如果 pitch 数据有更多列，可以根据需要调整
        pitch_value = pitch(matches, 4);
        match_indices = find(matches);
        % 如果有多个匹配，取第一个匹配的值
        pitch_result(i) = pitch_value(1);
    else
        % 如果找不到匹配，赋值为 0
        pitch_result(i) = 0;
    end
end

% 显示部分 pitch_result
numElementsToDisplay = min(10, length(pitch_result)); % 显示前10个元素
fprintf('前 %d 个 pitch_result 元素：\n', numElementsToDisplay);
for i = 1:numElementsToDisplay
    fprintf('pitch_result(%d) = %d\n', i, pitch_result(i));
end

% 如果需要，将 pitch_result 保存到文件
% save('pitch_result.mat', 'pitch_result');

% ------------------- 结束新增部分 -------------------
max3 = 3.0;  % heave 方向上的阈值为3米
max4 = 0.52;  % roll 30度
max5 = 0.52;  % pitch 30度

RAO_3 = heave_result;
RAO_4 = roll_result;
RAO_5 = pitch_result;

load('energy_block_0-255.mat')
energy_density =  squeeze(energy_density_half(element2, element1, :));
% 计算每个自由度的方差和协方差
var_x3 = sum(energy_density .* abs(RAO_3).^2) / (max3 * max3); % heave方差
var_x4 = sum(energy_density .* abs(RAO_4).^2) / (max4 * max4); % roll方差
var_x5 = sum(energy_density .* abs(RAO_5).^2) / (max5 * max5); % pitch方差
% 协方差计算公式：Cov(x_i, x_j) = sum(a_i * (RAO_i(ω) * RAO_j*(ω))^*)
cov_x3_x4 = sum(energy_density .* abs((RAO_3 .* conj(RAO_4)))) / (max3 * max4); % heave-roll协方差
cov_x3_x5 = sum(energy_density .* abs((RAO_3 .* conj(RAO_5)))) / (max3 * max5); % heave-pitch协方差
cov_x4_x5 = sum(energy_density .* abs((RAO_4 .* conj(RAO_5)))) / (max4 * max5); % roll-pitch协方差

% 定义最小阈值
min_threshold = 1e-9;

% 检查并调整每个变量，如果小于阈值则赋值为阈值
var_x3 = max(var_x3, min_threshold);
var_x4 = max(var_x4, min_threshold);
var_x5 = max(var_x5, min_threshold);
cov_x3_x4 = max(cov_x3_x4, min_threshold);
cov_x3_x5 = max(cov_x3_x5, min_threshold);
cov_x4_x5 = max(cov_x4_x5, min_threshold);
% 构造协方差矩阵
Sigma = [
    var_x3, cov_x3_x4, cov_x3_x5;
    cov_x3_x4, var_x4, cov_x4_x5;
    cov_x3_x5, cov_x4_x5, var_x5
];

% 输出协方差矩阵
disp('Covariance Matrix Sigma:');
disp(Sigma);

% 1. 使用协方差矩阵 Sigma
mu = [0, 0, 0];  % 0均值

% 2. 计算 Sigma 的 Cholesky 分解
L = chol(Sigma, 'lower');  % 'lower' 表示返回下三角矩阵

% 3. 生成标准正态分布样本
num_samples = 10000;
z = randn(num_samples, 3);  % 生成 10000 个标准正态分布样本

% 4. 通过 Cholesky 分解转换为具有协方差 Sigma 的样本
samples = z * L + mu;  % 变换样本，使其具有目标协方差

% 计算三个维度的最大值和最小值
max_val = max(abs(samples(:)));  % 取所有样本中绝对值的最大值
min_val = -max_val;  % 对称于原点

% 计算样本的协方差矩阵
cov_matrix = cov(samples);

% 计算协方差矩阵的特征值和特征向量
[EigVec, EigVal] = eig(cov_matrix);

% 特征值的平方根（标准差）
std_devs = sqrt(diag(EigVal));

% 计算95%置信椭球的尺度
confidence_interval = 1.96;  % 对应95%置信度

% 绘制 3D 散点图
figure;
scatter3(samples(:, 1), samples(:, 2), samples(:, 3), 5, 'filled', 'MarkerFaceAlpha', 0.3);

% 保证坐标轴范围相同
xlim([min_val, max_val]);
ylim([min_val, max_val]);
zlim([min_val, max_val]);

% 添加标签和标题
xlabel('Heave');
ylabel('Roll');
zlabel('Pitch');
title('3D Distribution with 95% Confidence Ellipsoid');
grid on;
axis equal;

% 绘制置信椭球
hold on;
% 计算95%置信椭球的参数并绘制
theta = linspace(0, 2*pi, 100);
phi = linspace(0, pi, 100);
[TH, PHI] = meshgrid(theta, phi);

% 将椭球的参数化表示转换为三维坐标
x_ellipsoid = confidence_interval * std_devs(1) * sin(PHI) .* cos(TH);
y_ellipsoid = confidence_interval * std_devs(2) * sin(PHI) .* sin(TH);
z_ellipsoid = confidence_interval * std_devs(3) * cos(PHI);

% 旋转椭球
ellipsoid_points = [x_ellipsoid(:), y_ellipsoid(:), z_ellipsoid(:)] * EigVec';
x_ellipsoid = reshape(ellipsoid_points(:, 1), size(x_ellipsoid));
y_ellipsoid = reshape(ellipsoid_points(:, 2), size(y_ellipsoid));
z_ellipsoid = reshape(ellipsoid_points(:, 3), size(z_ellipsoid));

% 绘制置信椭球
h = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 设置视角和显示比例
view(30, 30);

% 添加图例
hold on;
h_scatter = scatter3(samples(:, 1), samples(:, 2), samples(:, 3), 5, 'filled', 'MarkerFaceAlpha', 0.3);
h_ellipsoid = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 设置图例字体大小
%legend([h_scatter, h_ellipsoid], {'Samples', '95% Confidence Ellipsoid'}, 'Location', 'best', 'FontSize', 14);
legend([h_scatter, h_ellipsoid], {'Samples', '95% Confidence Ellipsoid'}, 'FontSize', 14, 'Position', [0.75, 0.75, 0.2, 0.2]);
% ------------------- 示例脚本: example_call_compute_cvar.m -------------------

% 清空工作区和命令窗口
%clear;
%clc;

% 1. 定义示例协方差矩阵 Sigma (3x3)
% Sigma = [
%     0.0000000005, 0.000000002, 0.0000000015;
%     0.000000002, 0.0000000004, 0.0000000008;
%     0.0000000015, 0.0000000008, 0.0000000003
% ];

% 2. 定义权重向量 w (3x1)
w = [1.0; 1.0; 1.0];  % 均匀权重

% 3. 定义置信水平 alpha
alpha = 0.95;  % 95% 置信水平

% 4. 指定输出文本文件的路径
%    确保该路径存在，或者MATLAB会报错
output_file = 'cvar_results_example.txt';  % 当前工作目录下

% 5. 调用 compute_cvar_and_save 函数
compute_cvar_and_save(Sigma, w, alpha, output_file);

% 6. 展示输出文件的内容
%    使用 MATLAB 的 fileread 函数读取并显示文本文件内容
fprintf('\n--- 输出文件内容 ---\n');
file_content = fileread(output_file);
disp(file_content);

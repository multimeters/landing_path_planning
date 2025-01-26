% ============================== 参数设置 ==============================

% 指定方向差异文件路径
directionFilePath = 'output_data/direction_differences_path_0.txt';

% 指定深度文件路径
depthFilePath = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/source/depth_a15_512x512.txt'; % 请将此路径替换为实际深度文件的路径

% 指定 omega.mat 文件路径
omegaMatPath = '\\wsl.localhost\Ubuntu/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/source/omega.mat';

% 指定 heave.mat、roll.mat、pitch.mat 文件路径
heaveMatPath = 'heave.mat'; % 请确保 heave.mat 的路径正确
rollMatPath = 'roll.mat';   % 请确保 roll.mat 的路径正确
pitchMatPath = 'pitch.mat'; % 请确保 pitch.mat 的路径正确

% 指定 energy_block 文件路径
energyBlockPath = 'energy_block_0-255.mat';

% 指定输出文件路径
outputFilePath = 'output_data/direction_differences_with_cvar.txt';

% 定义 CVaR 计算所需的其他参数
U_kmh = 10; % 速度 U，单位为 km/h
U = U_kmh * (1000/3600); % 将速度转换为 m/s
g = 9.81; % 重力加速度，单位为 m/s²
min_threshold = 1e-9; % 定义最小阈值

% 定义 compute_cvar_and_save 函数的权重向量和置信水平
w = [1.0; 1.0; 1.0];  % 均匀权重
alpha = 0.95;          % 95% 置信水平

% ============================== 检查文件是否存在 ==============================

% 检查方向差异文件是否存在
if ~exist(directionFilePath, 'file')
    error('文件 "%s" 未找到。请检查路径是否正确。', directionFilePath);
end

% 检查深度文件是否存在
if ~exist(depthFilePath, 'file')
    error('深度文件 "%s" 未找到。请检查路径是否正确。', depthFilePath);
end

% 检查 omega.mat 文件是否存在
if ~exist(omegaMatPath, 'file')
    error('omega.mat 文件 "%s" 未找到。请检查路径是否正确。', omegaMatPath);
end

% 检查 heave.mat 文件是否存在
if ~exist(heaveMatPath, 'file')
    error('heave.mat 文件 "%s" 未找到。请检查路径是否正确。', heaveMatPath);
end

% 检查 roll.mat 文件是否存在
if ~exist(rollMatPath, 'file')
    error('roll.mat 文件 "%s" 未找到。请检查路径是否正确。', rollMatPath);
end

% 检查 pitch.mat 文件是否存在
if ~exist(pitchMatPath, 'file')
    error('pitch.mat 文件 "%s" 未找到。请检查路径是否正确。', pitchMatPath);
end

% 检查 energy_block 文件是否存在
if ~exist(energyBlockPath, 'file')
    error('energy_block 文件 "%s" 未找到。请检查路径是否正确。', energyBlockPath);
end

% ============================== 加载数据 ==============================

% 加载深度数据
depthData = load(depthFilePath);
depthSize = size(depthData);
fprintf('深度文件 "%s" 的尺寸为：\n', depthFilePath);
disp(depthSize);

% 加载 omega.mat 文件
omegaData = load(omegaMatPath);
if isfield(omegaData, 'omega_1024')
    omega = omegaData.omega_1024;
    fprintf('成功加载 omega 数组，尺寸为：\n');
    disp(size(omega));
else
    error('omega.mat 中未找到变量 "omega_1024"。请检查变量名。');
end

% 加载 heave.mat 文件
heaveData = load(heaveMatPath);
if isfield(heaveData, 'heave')
    heave = heaveData.heave;
    fprintf('成功加载 heave 数据，尺寸为：\n');
    disp(size(heave));
else
    error('heave.mat 中未找到变量 "heave"。请检查变量名。');
end

% 加载 roll.mat 文件
rollData = load(rollMatPath);
if isfield(rollData, 'roll')
    roll = rollData.roll;
    fprintf('成功加载 roll 数据，尺寸为：\n');
    disp(size(roll));
else
    error('roll.mat 中未找到变量 "roll"。请检查变量名。');
end

% 加载 pitch.mat 文件
pitchData = load(pitchMatPath);
if isfield(pitchData, 'pitch')
    pitch = pitchData.pitch;
    fprintf('成功加载 pitch 数据，尺寸为：\n');
    disp(size(pitch));
else
    error('pitch.mat 中未找到变量 "pitch"。请检查变量名。');
end

% 加载 energy_block 文件
energyBlockData = load(energyBlockPath);
if isfield(energyBlockData, 'energy_density_half')
    energy_density_half = energyBlockData.energy_density_half;
    fprintf('成功加载 energy_density_half 数据，尺寸为：\n');
    disp(size(energy_density_half));
else
    error('energy_block 文件中未找到变量 "energy_density_half"。请检查变量名。');
end

% ============================== 读取方向差异文件 ==============================

% 读取方向差异文件
directionData = load(directionFilePath);
directionSize = size(directionData);
fprintf('文件 "%s" 的尺寸为：\n', directionFilePath);
disp(directionSize);

% 检查数据是否有至少5列
if directionSize(2) < 5
    error('方向差异数据的列数不足5列。');
end

% ============================== 初始化输出数据 ==============================

% 预分配输出数据矩阵，增加一列用于存储 CVaR
outputData = zeros(directionSize(1), 6); % 假设前5列数据，6列为 CVaR

% 将原始方向差异数据的前5列复制到输出数据
outputData(:, 1:5) = directionData(:, 1:5);

% ============================== 打开输出文件 ==============================

% 打开输出文件以写入
fid_out = fopen(outputFilePath, 'w');
if fid_out == -1
    error('无法打开输出文件 "%s" 进行写入。', outputFilePath);
end

% 如果需要写入表头，可以取消注释以下行并根据需要修改
% fprintf(fid_out, 'Column1 Column2 Column3 Column4 Column5 CVaR\n');

% ============================== 循环处理每一行数据 ==============================

% 遍历每一行数据
for row = 1:directionSize(1)
    fprintf('处理第 %d 行数据...\n', row);
    
    % 提取当前行的第1、2和5个数据
    element1 = directionData(row, 1); % 列索引
    element2 = directionData(row, 2); % 行索引
    element5 = directionData(row, 5); % μ
    
    % 验证 element1 和 element2 是否为整数
    if ~(isinteger(element1) || floor(element1) == element1) || ~(isinteger(element2) || floor(element2) == element2)
        warning('第 %d 行的 element1 或 element2 不是整数。跳过此行。', row);
        cvar = NaN;
        outputData(row, 6) = cvar;
        fprintf(fid_out, '%d %d %f %f %f %f\n', directionData(row, 1:5), cvar);
        continue;
    end
    
    % 转换为双精度以防止索引问题
    element1 = double(element1);
    element2 = double(element2);
    
    % 验证 element1 和 element2 是否在深度数据的范围内
    if element2 < 1 || element2 > depthSize(1) || element1 < 1 || element1 > depthSize(2)
        warning('第 %d 行的 element1 或 element2 超出了深度数据的范围。跳过此行。', row);
        cvar = NaN;
        outputData(row, 6) = cvar;
        fprintf(fid_out, '%d %d %f %f %f %f\n', directionData(row, 1:5), cvar);
        continue;
    end
    
    % 提取对应的深度 h
    h = depthData(element2, element1);
    
    % 检查 h 是否为正数
    if h <= 0
        warning('第 %d 行的深度 h 必须为正数。当前 h = %f。跳过此行。', row, h);
        cvar = NaN;
        outputData(row, 6) = cvar;
        fprintf(fid_out, '%d %d %f %f %f %f\n', directionData(row, 1:5), cvar);
        continue;
    end
    
    % 假设 μ (element5) 是以度为单位，如果是弧度请跳过此步骤
    mu_degrees = element5;
    mu = deg2rad(mu_degrees); % 将 μ 转换为弧度
    
    % 计算表达式 1 - (U * cos(mu)) / sqrt(g * h)
    term = (U * cos(mu)) / sqrt(g * h);
    result = 1 - term;
    
    % 计算修正频率 omega_corrected
    omega_corrected = omega * result;
    
    % 处理 omega_corrected:
    % 1. 每个元素除以0.1
    % 2. 四舍五入取整
    % 3. 如果结果为 Inf 或 -Inf，设为0
    omega_processed = omega_corrected / 0.1;
    omega_processed = round(omega_processed);
    inf_indices = isinf(omega_corrected);
    omega_processed(inf_indices) = 0;
    
    % 处理 mu
    mu_processed = mu / 0.5;
    mu_processed = round(mu_processed);
    
    % 提取 heave_result
    heave_result = zeros(length(omega_processed), 1);
    for i = 1:length(omega_processed)
        current_omega = omega_processed(i);
        current_mu = mu_processed;
        matches = (heave(:,1) == current_omega) & (heave(:,2) == current_mu);
        if any(matches)
            heave_value = heave(matches, 4);
            heave_result(i) = heave_value(1);
        else
            heave_result(i) = 0;
        end
    end
    
    % 提取 roll_result
    roll_result = zeros(length(omega_processed), 1);
    for i = 1:length(omega_processed)
        current_omega = omega_processed(i);
        current_mu = mu_processed;
        matches = (roll(:,1) == current_omega) & (roll(:,2) == current_mu);
        if any(matches)
            roll_value = roll(matches, 4);
            roll_result(i) = roll_value(1);
        else
            roll_result(i) = 0;
        end
    end
    
    % 提取 pitch_result
    pitch_result = zeros(length(omega_processed), 1);
    for i = 1:length(omega_processed)
        current_omega = omega_processed(i);
        current_mu = mu_processed;
        matches = (pitch(:,1) == current_omega) & (pitch(:,2) == current_mu);
        if any(matches)
            pitch_value = pitch(matches, 4);
            pitch_result(i) = pitch_value(1);
        else
            pitch_result(i) = 0;
        end
    end
    
    % 定义阈值
    max3 = 3.0;  % heave 方向上的阈值为3米
    max4 = 0.52; % roll 30度
    max5 = 0.52; % pitch 30度
    
    RAO_3 = heave_result;
    RAO_4 = roll_result;
    RAO_5 = pitch_result;
    
    % 提取 energy_density
    energy_density = squeeze(energy_density_half(element2, element1, :));
    
    % 检查 energy_density 的维度是否匹配
    if length(energy_density) ~= length(RAO_3)
        warning('第 %d 行的 energy_density 长度与 RAO_3 不匹配。跳过此行。', row);
        cvar = NaN;
        outputData(row, 6) = cvar;
        fprintf(fid_out, '%d %d %f %f %f %f\n', directionData(row, 1:5), cvar);
        continue;
    end
    
    % 计算方差和协方差
    var_x3 = sum(energy_density .* abs(RAO_3).^2) / (max3^2); % heave方差
    var_x4 = sum(energy_density .* abs(RAO_4).^2) / (max4^2); % roll方差
    var_x5 = sum(energy_density .* abs(RAO_5).^2) / (max5^2); % pitch方差
    cov_x3_x4 = sum(energy_density .* (RAO_3 .* RAO_4)) / (max3 * max4); % heave-roll协方差
    cov_x3_x5 = sum(energy_density .* (RAO_3 .* RAO_5)) / (max3 * max5); % heave-pitch协方差
    cov_x4_x5 = sum(energy_density .* (RAO_4 .* RAO_5)) / (max4 * max5); % roll-pitch协方差
    
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
    
    % 使用预先编写好的 compute_cvar_and_save 函数计算 CVaR
    % 假设 compute_cvar_and_save(Sigma, w, alpha, output_file) 会返回 CVaR 值
    try
        cvar = compute_cvar_and_save(Sigma, w, alpha, outputFilePath); % 修改为直接返回 CVaR
    catch ME
        warning('计算第 %d 行的 CVaR 时出错：%s', row, ME.message);
        cvar = NaN;
    end
    
    % 检查 cvar 是否有效
    if isempty(cvar) || ~isnumeric(cvar)
        warning('第 %d 行的 CVaR 为空或非数值。将 CVaR 设为 NaN。', row);
        cvar = NaN;
    end
    
    % 将 CVaR 添加到输出数据中
    outputData(row, 6) = cvar;
    
    % 将当前行数据和 CVaR 写入输出文件
    % 格式: 前5列数据保持不变，最后一列为 CVaR
    fprintf(fid_out, '%d %d %f %f %f %f\n', directionData(row, 1:5), cvar);
end

% ============================== 关闭输出文件 ==============================

fclose(fid_out);
fprintf('所有数据已处理完毕，并保存到 "%s"。\n', outputFilePath);

% ============================== 可选：保存输出数据为 MAT 文件 ==============================
% save('direction_differences_with_cvar.mat', 'outputData');

% ============================== 结束 ==============================

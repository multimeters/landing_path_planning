% ============================== 参数设置 ==============================

% 指定路径文件夹
inputFolder = 'output_data_path_e_01_lam_1001_sam_1000000/';

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

% 定义 CVaR 计算所需的其他参数
U_kmh = 10; % 速度 U，单位为 km/h
U = U_kmh * (1000/3600); % 将速度转换为 m/s
g = 9.81; % 重力加速度，单位为 m/s²
min_threshold = 1e-9; % 定义最小阈值

% 定义权重向量和置信水平
w = [1.0; 1.0; 1.0];  % 均匀权重
alpha = 0.95;          % 95% 置信水平

% ============================== 检查必要文件是否存在 ==============================

% 列出所有需要检查的文件及其描述
requiredFiles = {
    depthFilePath, '深度文件';
    omegaMatPath, 'omega.mat 文件';
    heaveMatPath, 'heave.mat 文件';
    rollMatPath, 'roll.mat 文件';
    pitchMatPath, 'pitch.mat 文件';
    energyBlockPath, 'energy_block 文件'
};

% 循环检查每个文件是否存在
for i = 1:size(requiredFiles, 1)
    filePath = requiredFiles{i, 1};
    fileDesc = requiredFiles{i, 2};
    if ~exist(filePath, 'file')
        error('文件 "%s" (%s) 未找到。请检查路径是否正确。', filePath, fileDesc);
    end
end

% ============================== 加载共享数据 ==============================

fprintf('加载深度数据...\n');
depthData = load(depthFilePath);
depthSize = size(depthData);
fprintf('深度文件 "%s" 的尺寸为：%d 行 x %d 列\n', depthFilePath, depthSize(1), depthSize(2));

fprintf('加载 omega.mat 文件...\n');
omegaData = load(omegaMatPath);
if isfield(omegaData, 'omega_1024')
    omega = omegaData.omega_1024;
    fprintf('成功加载 omega 数组，尺寸为：%d\n', numel(omega));
else
    error('omega.mat 中未找到变量 "omega_1024"。请检查变量名。');
end

fprintf('加载 heave.mat 文件...\n');
heaveData = load(heaveMatPath);
if isfield(heaveData, 'heave')
    heave = heaveData.heave;
    fprintf('成功加载 heave 数据，尺寸为：%d 行 x %d 列\n', size(heave, 1), size(heave, 2));
else
    error('heave.mat 中未找到变量 "heave"。请检查变量名。');
end

fprintf('加载 roll.mat 文件...\n');
rollData = load(rollMatPath);
if isfield(rollData, 'roll')
    roll = rollData.roll;
    fprintf('成功加载 roll 数据，尺寸为：%d 行 x %d 列\n', size(roll, 1), size(roll, 2));
else
    error('roll.mat 中未找到变量 "roll"。请检查变量名。');
end

fprintf('加载 pitch.mat 文件...\n');
pitchData = load(pitchMatPath);
if isfield(pitchData, 'pitch')
    pitch = pitchData.pitch;
    fprintf('成功加载 pitch 数据，尺寸为：%d 行 x %d 列\n', size(pitch, 1), size(pitch, 2));
else
    error('pitch.mat 中未找到变量 "pitch"。请检查变量名。');
end

fprintf('加载 energy_block 文件...\n');
energyBlockData = load(energyBlockPath);
if isfield(energyBlockData, 'energy_density_half')
    energy_density_half = energyBlockData.energy_density_half;
    fprintf('成功加载 energy_density_half 数据，尺寸为：%d x %d x %d\n', size(energy_density_half,1), size(energy_density_half,2), size(energy_density_half,3));
else
    error('energy_block 文件中未找到变量 "energy_density_half"。请检查变量名。');
end

% ============================== 创建查找表 ==============================

fprintf('创建 heave、roll 和 pitch 的查找表...\n');

% 创建查找表函数
createLookupMap = @(data) containers.Map(...
    strcat(string(data(:,1)), '_', string(data(:,2))), ...
    data(:,4));

% 创建 heave、roll 和 pitch 的查找表
heaveMap = createLookupMap(heave);
rollMap = createLookupMap(roll);
pitchMap = createLookupMap(pitch);

fprintf('查找表创建完成。\n');

% ============================== 计算常数 ==============================

fprintf('计算 z_alpha 和 pdf_z_alpha...\n');
z_alpha = norminv(alpha);
pdf_z_alpha = normpdf(z_alpha);
fprintf('z_alpha = %.6f, pdf_z_alpha = %.6f\n', z_alpha, pdf_z_alpha);

% ============================== 启动并行池 ==============================

fprintf('启动并行池（如果尚未启动）...\n');
pool = gcp('nocreate'); % 获取当前并行池
if isempty(pool)
    parpool; % 启动默认并行池
    fprintf('并行池已启动。\n');
else
    fprintf('使用已存在的并行池。\n');
end

% ============================== 批量处理所有路径文件 ==============================

% 获取所有方向差异路径文件
pathFiles = dir(fullfile(inputFolder, 'direction_differences_path_*.txt'));

% 检查是否找到任何路径文件
if isempty(pathFiles)
    error('在文件夹 "%s" 中未找到任何匹配的路径文件。', inputFolder);
end

% 循环处理每个路径文件
for pf = 1:length(pathFiles)
    % 获取当前文件的信息
    currentFile = pathFiles(pf);
    inputFilePath = fullfile(inputFolder, currentFile.name);
    
    % 提取路径编号以匹配输出文件名
    tokens = regexp(currentFile.name, 'direction_differences_path_(\d+)\.txt', 'tokens');
    if isempty(tokens)
        warning('文件 "%s" 不符合预期的命名格式。跳过此文件。', currentFile.name);
        continue;
    end
    pathNumber = tokens{1}{1};
    
    % 生成对应的输出文件名
    outputFileName = sprintf('direction_differences_with_cvar_path_%s.txt', pathNumber);
    outputFilePath = fullfile(inputFolder, outputFileName);
    
    fprintf('处理文件: %s -> %s\n', currentFile.name, outputFileName);
    
    try
        % 调用处理函数
        processPathFile(inputFilePath, outputFilePath, depthData, depthSize, ...
            omega, heaveMap, rollMap, pitchMap, energy_density_half, ...
            U, g, min_threshold, w, alpha, z_alpha, pdf_z_alpha);
        
        fprintf('文件 "%s" 处理完成，输出保存为 "%s"\n', currentFile.name, outputFileName);
    catch ME
        warning('处理文件 "%s" 时发生错误: %s', currentFile.name, ME.message);
    end
end

fprintf('所有路径文件处理完毕。\n');

% ============================== 结束 ==============================

% ============================== 处理函数定义 ==============================

function processPathFile(directionFilePath, outputFilePath, depthData, depthSize, ...
    omega, heaveMap, rollMap, pitchMap, energy_density_half, ...
    U, g, min_threshold, w, alpha, z_alpha, pdf_z_alpha)

    % ============================== 读取方向差异文件 ==============================
    
    fprintf('读取方向差异文件: %s\n', directionFilePath);
    directionData = load(directionFilePath);
    directionSize = size(directionData);
    fprintf('文件 "%s" 的尺寸为：%d 行 x %d 列\n', directionFilePath, directionSize(1), directionSize(2));
    
    % 检查数据是否有至少5列
    if directionSize(2) < 5
        error('方向差异数据的列数不足5列。');
    end
    
    % ============================== 初始化输出数据 ==============================
    
    fprintf('初始化输出数据...\n');
    % 预分配输出数据矩阵，增加一列用于存储 CVaR
    outputData = zeros(directionSize(1), 6); % 前5列数据，加1列 CVaR
    outputData(:, 1:5) = directionData(:, 1:5);
    
    % ============================== 并行处理数据 ==============================
    
    fprintf('开始并行处理数据...\n');
    
    % 转换 directionData 的第1和第2列为双精度以确保索引正确
    directionIndices1 = double(directionData(:,1));
    directionIndices2 = double(directionData(:,2));
    
    % 提取 mu_degrees 并转换为弧度
    mu_degrees = directionData(:,5);
    mu = deg2rad(mu_degrees);
    
    % 提取对应的深度 h
    h = depthData(sub2ind(depthSize, directionIndices2, directionIndices1));
    
    % 预计算 term 和 result
    term = (U * cos(mu)) ./ sqrt(g * h);
    result = 1 - term;
    
    % 计算 omega_corrected 并处理
    omega_corrected = omega(:) * result(:)'; % 形成矩阵，行对应 omega，列对应 directionData
    omega_processed = round(omega_corrected / 0.1);
    omega_processed(isinf(omega_corrected)) = 0;
    
    % 处理 mu
    mu_processed = round(mu / 0.5);
    
    % 获取 energy_density for all directions
    % Assuming energy_density_half is of size [rows, cols, omega_length]
    % We will reshape it for easier access in parallel
    energy_density_all = reshape(energy_density_half, [depthSize(1), depthSize(2), size(energy_density_half,3)]);
    
    % 定义阈值
    max3 = 3.0;  % heave 方向上的阈值为3米
    max4 = 0.52; % roll 30度
    max5 = 0.52; % pitch 30度
    
    % 预分配 CVaR 列
    cvar_all = NaN(directionSize(1), 1);
    
    % 使用 parfor 并行循环处理每一行数据
    parfor row = 1:directionSize(1)
        % 提取当前行的索引和 mu
        idx1 = directionIndices1(row);
        idx2 = directionIndices2(row);
        mu_deg = mu_degrees(row);
        mu_rad = mu(row);
        
        % 获取 h，已预计算
        current_h = h(row);
        
        % 检查 h 是否为正数
        if current_h <= 0
            warning('第 %d 行的深度 h 必须为正数。当前 h = %f。跳过此行。', row, current_h);
            cvar_all(row) = NaN;
            continue;
        end
        
        % 计算 term 和 result，已预计算
        current_result = result(row);
        
        % 计算 omega_corrected
        current_omega_corrected = omega * current_result;
        
        % 处理 omega_corrected:
        % 1. 每个元素除以0.1
        % 2. 四舍五入取整
        % 3. 如果结果为 Inf 或 -Inf，设为0
        current_omega_processed = round(current_omega_corrected / 0.1);
        current_omega_processed(isinf(current_omega_corrected)) = 0;
        
        % 处理 mu
        current_mu_processed = round(mu_rad / 0.5);
        
        % 构建键以查询查找表
        keys = strcat(string(current_omega_processed), '_', string(current_mu_processed));
        
        % 提取 heave_result、roll_result、pitch_result
        heave_result = zeros(length(current_omega_processed),1);
        roll_result = zeros(length(current_omega_processed),1);
        pitch_result = zeros(length(current_omega_processed),1);
        
        for i = 1:length(current_omega_processed)
            key = keys(i);
            if heaveMap.isKey(key)
                heave_result(i) = heaveMap(key);
            else
                heave_result(i) = 0;
            end
            
            if rollMap.isKey(key)
                roll_result(i) = rollMap(key);
            else
                roll_result(i) = 0;
            end
            
            if pitchMap.isKey(key)
                pitch_result(i) = pitchMap(key);
            else
                pitch_result(i) = 0;
            end
        end
        
        % 提取 energy_density
        current_energy_density = squeeze(energy_density_all(idx2, idx1, :));
        current_energy_density = current_energy_density(:);
        
        % 检查 energy_density 的长度是否与 RAO_* 一致
        if length(current_energy_density) ~= length(heave_result)
            warning('第 %d 行的 energy_density 长度与 RAO_* 不匹配。跳过此行。', row);
            cvar_all(row) = NaN;
            continue;
        end
        
        % 计算方差和协方差
        var_x3 = sum(current_energy_density .* (abs(heave_result).^2)) / (max3^2);
        var_x4 = sum(current_energy_density .* (abs(roll_result).^2)) / (max4^2);
        var_x5 = sum(current_energy_density .* (abs(pitch_result).^2)) / (max5^2);
        cov_x3_x4 = sum(current_energy_density .* (heave_result .* roll_result)) / (max3 * max4);
        cov_x3_x5 = sum(current_energy_density .* (heave_result .* pitch_result)) / (max3 * max5);
        cov_x4_x5 = sum(current_energy_density .* (roll_result .* pitch_result)) / (max4 * max5);
        
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
        
        % 检查 Sigma 是否为正定矩阵
        eigenvalues = eig(Sigma);
        if any(eigenvalues <= 0)
            warning('第 %d 行协方差矩阵不是正定的。跳过此行。', row);
            cvar_all(row) = NaN;
            continue;
        end
        
        % 计算 sigma_L
        sigma_L = w' * Sigma * w;
        
        % 计算 CVaR
        current_cvar = (pdf_z_alpha / (1 - alpha)) * sqrt(sigma_L);
        
        % 检查 cvar 是否有效
        if isempty(current_cvar) || ~isnumeric(current_cvar)
            warning('第 %d 行的 CVaR 为空或非数值。将 CVaR 设为 NaN。', row);
            current_cvar = NaN;
        end
        
        % 存储 CVaR
        cvar_all(row) = current_cvar;
    end
    
    fprintf('并行处理完成。\n');
    
    % ============================== 写入输出文件 ==============================
    
    fprintf('写入输出文件: %s\n', outputFilePath);
    
    % 将 CVaR 添加到输出数据中
    outputData(:,6) = cvar_all;
    
    % 创建格式化字符串
    % 前两列为整数，后四列为浮点数
    formatSpec = '%d %d %.6f %.6f %.6f %.6f\n';
    
    % 打开输出文件以写入
    fid_out = fopen(outputFilePath, 'w');
    if fid_out == -1
        error('无法打开输出文件 "%s" 进行写入。', outputFilePath);
    end
    
    % 写入所有数据
    % 使用矢量化的 fprintf
    fprintf(fid_out, formatSpec, outputData.');
    
    % 关闭输出文件
    fclose(fid_out);
    fprintf('输出文件 "%s" 已成功写入。\n', outputFilePath);
    
end

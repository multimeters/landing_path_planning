% 定义文件路径
uInputDir = 'pixel_data/u'; % u 文件目录
vInputDir = 'pixel_data/v'; % v 文件目录
outputDir = 'pixel_data/a'; % 输出文件目录

% 初始化一些参数
numFiles = 1000; % 文件数量
blockSize = 256; % 每块的行数
numBlocks = 4096 / blockSize; % 块的数量

% 开始计时
tic;
disp('程序开始执行...');

% 遍历所有块
for blockIdx = 1:numBlocks
    blockIdx
    % 初始化存储角度数据的矩阵
    angles = zeros(blockSize, 2048, numFiles);
    
    % 遍历所有文件
    for fileIdx = 0:(numFiles-1)
        fileIdx
        % 构建文件名
        uFileName = sprintf('u_%05d.mat', fileIdx);
        vFileName = sprintf('v_%05d.mat', fileIdx);
        uFilePath = fullfile(uInputDir, uFileName);
        vFilePath = fullfile(vInputDir, vFileName);
        
        % 加载 u 和 v 文件并提取数据
        loadedUData = load(uFilePath, 'data');  % 仅加载所需变量
        loadedVData = load(vFilePath, 'data');
        u = loadedUData.data;
        v = loadedVData.data;
        
        % 提取块
        startRow = (blockIdx - 1) * blockSize + 1;
        endRow = blockIdx * blockSize;
        uBlock = u(startRow:endRow, :);
        vBlock = v(startRow:endRow, :);
        
        % 计算每个点的角度（以弧度为单位）
        a(:, :, fileIdx+1) = atan2(vBlock, uBlock);
        
        % 清除加载的数据以释放内存
        clear u v loadedUData loadedVData uBlock vBlock;

        % 输出当前进度信息
        if mod(fileIdx + 1, 100) == 0 || fileIdx == numFiles - 1
            fprintf('已处理文件：%d/%d，正在处理块：%d/%d\n', fileIdx + 1, numFiles, blockIdx, numBlocks);
        end
    end

    % 构建输出文件名，格式为 angle_startRow_endRow.mat
    outputFileName = sprintf('a_%d_%d.mat', startRow, endRow);
    outputPath = fullfile(outputDir, outputFileName);
    
    % 保存每块的角度数据
    save(outputPath, 'a', '-v7.3');

    % 清除角度数据以释放内存
    clear a;
    
    % 输出保存文件的信息
    fprintf('块 %d (%d-%d) 已保存为 %s\n', blockIdx, startRow, endRow, outputFileName);
end

% 显示总执行时间
elapsedTime = toc;
fprintf('程序执行完毕，总耗时：%.2f 秒\n', elapsedTime);

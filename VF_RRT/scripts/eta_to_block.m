% 定义文件路径
inputDir = ''; % 输入文件目录
outputDir = ''; % 输出文件目录

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
    % 初始化存储块数据的矩阵
    data = zeros(blockSize, 2048, numFiles);
    
    % 遍历所有文件
    for fileIdx = 0:(numFiles-1)
        fileIdx
        % 构建文件名
        fileName = sprintf('eta_%05d.mat', fileIdx);
        filePath = fullfile(inputDir, fileName);
        
        % 加载文件并提取数据
        loadedData = load(filePath, 'data');  % 仅加载所需变量
        eta = loadedData.data;
        
        % 提取块
        startRow = (blockIdx - 1) * blockSize + 1;
        endRow = blockIdx * blockSize;
        data(:, :, fileIdx+1) = eta(startRow:endRow, :);
        
        % 清除加载的数据以释放内存
        clear eta loadedData;

        % 输出当前进度信息
        if mod(fileIdx + 1, 100) == 0 || fileIdx == numFiles - 1
            fprintf('已处理文件：%d/%d，正在处理块：%d/%d\n', fileIdx + 1, numFiles, blockIdx, numBlocks);
        end
    end

    % 构建输出文件名，格式为 eta_startRow_endRow.mat
    outputFileName = sprintf('eta_%d_%d.mat', startRow, endRow);
    outputPath = fullfile(outputDir, outputFileName);
    
    % 保存每块的数据
    % 保存每块的数据
    save(outputPath, 'data', '-v7.3');


    % 清除块数据以释放内存
    clear dataBlock;
    
    % 输出保存文件的信息
    fprintf('块 %d (%d-%d) 已保存为 %s\n', blockIdx, startRow, endRow, outputFileName);
end

% 显示总执行时间
elapsedTime = toc;
fprintf('程序执行完毕，总耗时：%.2f 秒\n', elapsedTime);

% 定义测试程序参数
inputDirBefore = ''; % 转换前的文件目录
inputDirAfter = '';  % 转换后的文件目录
numFiles = 1000; % 文件数量
blockSize = 256; % 每块的行数
numBlocks = 4096 / blockSize; % 块的数量
numSamples = 100; % 采样数

% 随机生成采样帧数和像素位置
sampleFrames = randi([0 numFiles-1], numSamples, 1); % 随机生成帧数
samplePositions = [randi([1 4096], numSamples, 1), randi([1 2048], numSamples, 1)]; % 随机生成像素位置 [row, col]

% 进行一致性测试
for i = 1:numSamples
    i
    % 获取采样信息
    frameIdx = sampleFrames(i);
    rowIdx = samplePositions(i, 1);
    colIdx = samplePositions(i, 2);
    
    % 计算块索引
    blockIdx = ceil(rowIdx / blockSize);
    startRow = (blockIdx - 1) * blockSize + 1;
    endRow = blockIdx * blockSize;
    
    % 在块中的行索引
    blockRowIdx = rowIdx - startRow + 1;
    
    % 加载转换前的数据
    fileNameBefore = sprintf('eta_%05d.mat', frameIdx);
    filePathBefore = fullfile(inputDirBefore, fileNameBefore);
    loadedDataBefore = load(filePathBefore, 'data');
    etaBefore = loadedDataBefore.data(rowIdx, colIdx);
    
    % 加载转换后的数据
    fileNameAfter = sprintf('eta_%d_%d.mat', startRow, endRow);
    filePathAfter = fullfile(inputDirAfter, fileNameAfter);
    loadedDataAfter = load(filePathAfter, 'data');
    etaAfter = loadedDataAfter.data(blockRowIdx, colIdx, frameIdx + 1);
    
    % 对比数据是否一致
    if etaBefore ~= etaAfter
        fprintf('不一致发现：帧数 %d，位置 (%d, %d)，转换前数据 %.6f，转换后数据 %.6f\n', frameIdx, rowIdx, colIdx, etaBefore, etaAfter);
    else
        fprintf('一致：帧数 %d，位置 (%d, %d)\n', frameIdx, rowIdx, colIdx);
    end
    
    % 清除加载的数据以释放内存
    clear etaBefore etaAfter loadedDataBefore loadedDataAfter;
end

fprintf('测试完毕。\n');

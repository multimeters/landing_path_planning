% 定义测试程序参数
uInputDir = 'pixel_data/u'; % u 文件目录
vInputDir = 'pixel_data/v'; % v 文件目录
aInputDir = 'pixel_data/a'; % a 文件目录
numFiles = 1000; % 文件数量
blockSize = 256; % 每块的行数
numBlocks = 4096 / blockSize; % 块的数量
numSamples = 100; % 采样数

% 随机生成采样帧数和像素位置
sampleFrames = randi([0 numFiles-1], numSamples, 1); % 随机生成帧数
samplePositions = [randi([1 4096], numSamples, 1), randi([1 2048], numSamples, 1)]; % 随机生成像素位置 [row, col]

% 进行一致性测试
for i = 1:numSamples
    % 获取采样信息
    i
    frameIdx = sampleFrames(i);
    rowIdx = samplePositions(i, 1);
    colIdx = samplePositions(i, 2);
    
    % 计算块索引
    blockIdx = ceil(rowIdx / blockSize);
    startRow = (blockIdx - 1) * blockSize + 1;
    endRow = blockIdx * blockSize;
    
    % 在块中的行索引
    blockRowIdx = rowIdx - startRow + 1;
    
    % 加载原始 u 和 v 数据
    uFileName = sprintf('u_%05d.mat', frameIdx);
    vFileName = sprintf('v_%05d.mat', frameIdx);
    uFilePath = fullfile(uInputDir, uFileName);
    vFilePath = fullfile(vInputDir, vFileName);
    
    loadedUData = load(uFilePath, 'data');
    loadedVData = load(vFilePath, 'data');
    
    uValue = loadedUData.data(rowIdx, colIdx);
    vValue = loadedVData.data(rowIdx, colIdx);
    
    % 计算角度
    calculatedAngle = atan2(vValue, uValue);
    
    % 加载转换后的角度数据
    aFileName = sprintf('a_%d_%d.mat', startRow, endRow);
    aFilePath = fullfile(aInputDir, aFileName);
    loadedAData = load(aFilePath, 'a');
    storedAngle = loadedAData.a(blockRowIdx, colIdx, frameIdx + 1);
    
    % 对比数据是否一致
    if abs(calculatedAngle - storedAngle) > 1e-6
        fprintf('不一致发现：帧数 %d，位置 (%d, %d)，计算角度 %.6f，存储角度 %.6f\n', frameIdx, rowIdx, colIdx, calculatedAngle, storedAngle);
    else
        fprintf('一致：帧数 %d，位置 (%d, %d)\n', frameIdx, rowIdx, colIdx);
    end
    
    % 清除加载的数据以释放内存
    clear loadedUData loadedVData loadedAData uValue vValue calculatedAngle storedAngle;
end

fprintf('测试完毕。\n');

% 加载数据
load('encountered_frequency_coefficient_4096_2048.mat', 'encountered_frequency_coefficient');
load('heave_motionRAO_at_0.mat', 'motionRAO_w', 'motionRAO_amp_at_0');

% 设置参数
N = 1000;
fs = 10; % 设置采样频率
omega = linspace(0, 2*pi*fs, N);

% 指定像素坐标
pixel_row = 1000;
pixel_col = 1592;

% 计算对应的文件和行列索引
file_idx = ceil(pixel_row / 256);
row_in_file = mod(pixel_row - 1, 256) + 1;

% 加载对应的.mat文件
file_name = sprintf('pixel_data/eta/eta_%d_%d.mat', (file_idx-1)*256+1, file_idx*256);
data = load(file_name);

% 获取对应像素的时间序列数据
pixel_data = squeeze(data.data(row_in_file, pixel_col, :));

% 计算FFT
N = length(pixel_data);
fft_result = fft(pixel_data);
fft_result = fft_result(1:floor(N/2)+1);

% 计算能量密度函数
psd = (1/(2*pi*N)) * abs(fft_result).^2;
psd(2:end-1) = 2*psd(2:end-1);

% 计算频率向量
frequencies = (fs*(0:(N/2)))/N;

% 查找能量密度大于0.1的频率
significant_indices = find(psd > 0.1);
significant_frequencies = frequencies(significant_indices);

% 获取encountered_frequency_coefficient参数
coeff = encountered_frequency_coefficient(pixel_row, pixel_col);

% 检查coeff是否为NaN
if isnan(coeff)
    total_sum = NaN;
else
    % 修正非零频率
    adjusted_frequencies = significant_frequencies * coeff;

    % 查找motionRAO_w中最接近adjusted_frequencies的索引
    closest_indices = arrayfun(@(x) find(abs(motionRAO_w-x) == min(abs(motionRAO_w-x)), 1, 'first'), adjusted_frequencies);

    % 获取motionRAO_amp_at_0中对应的幅度值
    resulting_amplitudes = motionRAO_amp_at_0(closest_indices)';

    % 计算幅度的平方与对应能量密度的乘积
    amplitude_squared = resulting_amplitudes .^ 2;
    significant_psd_values = psd(significant_indices);
    product = amplitude_squared .* significant_psd_values;

    % 计算总和
    total_sum = sum(product);
end

% 绘制能量密度函数图
figure;
plot(frequencies, psd);
xlabel('Frequency (Hz)');
ylabel('Energy Density');
title(sprintf('Energy Density Function at pixel (%d, %d)', pixel_row, pixel_col));
grid on;

% 输出结果
disp('Total sum of products:');
disp(total_sum);

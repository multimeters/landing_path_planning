% 加载数据
load('encountered_frequency_coefficient_4096_2048.mat', 'encountered_frequency_coefficient');
load('heave_motionRAO_at_0.mat', 'motionRAO_w', 'motionRAO_amp_at_0');

% 设置参数
N = 1000;
fs = 10; % 设置采样频率

% 初始化总和矩阵
total_sum_matrix = zeros(4096, 2048);

% 关闭现有并行池（如果存在）
% if ~isempty(gcp('nocreate'))
%     delete(gcp('nocreate'));
% end
% 
% % 创建具有2个工作器的并行池
% parpool('local', 2);
tic; % 开始计时

% 遍历所有16个文件
for file_idx = 1:16
    % 加载对应的.mat文件
    file_name = sprintf('pixel_data/eta/eta_%d_%d.mat', (file_idx-1)*256+1, file_idx*256);
    data = load(file_name, 'data');
    
    % 使用for处理每一行
    for row = 1:256
        % 当前全局行号
        global_row = (file_idx-1)*256 + row
        
        for col = 1:2048
            % 获取encountered_frequency_coefficient参数
            coeff = encountered_frequency_coefficient(global_row, col);
            
            % 如果coeff为NaN，则直接将结果设为NaN
            if isnan(coeff)
                total_sum_matrix(global_row, col) = NaN;
                continue; % 跳过本次循环
            end

            % 获取对应像素的时间序列数据
            pixel_data = squeeze(data.data(row, col, :));

            % 计算FFT
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

            % 修正非零频率
            adjusted_frequencies = significant_frequencies * coeff;

            % 查找motionRAO_w中最接近adjusted_frequencies的索引
            closest_indices_cell = arrayfun(@(x) find(abs(motionRAO_w-x) == min(abs(motionRAO_w-x)), 1, 'first'), adjusted_frequencies, 'UniformOutput', false);
            closest_indices = cellfun(@(c) c(1), closest_indices_cell);

            % 获取motionRAO_amp_at_0中对应的幅度值
            resulting_amplitudes = motionRAO_amp_at_0(closest_indices)';

            % 计算幅度的平方与对应能量密度的乘积
            amplitude_squared = resulting_amplitudes .^ 2;
            significant_psd_values = psd(significant_indices);
            product = amplitude_squared .* significant_psd_values;

            % 计算总和并累加到total_sum_matrix中
            total_sum_matrix(global_row, col) = sum(product);
        end
    end
end

toc % 结束计时并保存该行的执行时间
% 显示处理完成信息
disp('Processing complete. Total sums calculated for each pixel.');

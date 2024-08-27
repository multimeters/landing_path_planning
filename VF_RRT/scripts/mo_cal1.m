% 设置参数
N = 1000;
fs = 10; % 请根据您的具体情况设置采样频率
omega = linspace(0, 2*pi*fs, N);

% 取 omega 小于 5 的部分
omega_half = omega(omega < 5);

% 加载数据
load('encountered_frequency_coefficient_4096_2048.mat', 'encountered_frequency_coefficient');
load('heave_motionRAO_at_0.mat', 'motionRAO_w', 'motionRAO_amp_at_0');

% 分块处理参数
block_size = 256;
num_blocks = size(encountered_frequency_coefficient, 1) / block_size;

for block = 1:num_blocks
    % 确定当前块的范围
    start_row = (block - 1) * block_size + 1;
    end_row = block * block_size;
    
    % 提取当前块的数据
    current_block = encountered_frequency_coefficient(start_row:end_row, :);
    
    % 初始化当前块的结果
    current_result = zeros(block_size, 2048, length(omega_half));
    
    % 对于每个像素，处理当前块的数据
    for i = 1:block_size
        i
        for j = 1:size(encountered_frequency_coefficient, 2)
            j
            if(~isnan(current_block(i, j)))
                encountered_omega_half = current_block(i, j) * omega_half;
                indices = arrayfun(@(freq) find(abs(motionRAO_w - freq) == min(abs(motionRAO_w - freq)), 1), encountered_omega_half, 'UniformOutput', false);
                indices = cell2mat(indices);
                % 检查索引是否为空
                if ~isempty(indices)
                    current_result(i, j, :) = motionRAO_amp_at_0(indices).^2;
                else
                    current_result(i, j, :) = NaN; % 或者设置为某个默认值
                end
            else
                current_result(i, j, :) = NaN; % 或者设置为某个默认值
            end
        end
    end
    
    % 保存当前块的结果
    save_filename = sprintf('processed_data_%d_%d.mat', start_row, end_row);
    save(save_filename, 'current_result');
end

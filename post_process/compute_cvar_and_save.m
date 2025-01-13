function Cvar=compute_cvar_and_save(Sigma, w, alpha, output_file)
    % compute_cvar_and_save 计算条件风险价值（CVaR）并将结果保存到文本文件
    %
    % 语法:
    %   compute_cvar_and_save(Sigma, w, alpha, output_file)
    %
    % 输入参数:
    %   Sigma       - 3x3 协方差矩阵
    %   w           - 3x1 权重向量
    %   alpha       - 置信水平 (0 < alpha < 1)
    %   output_file - 字符串，结果保存的文本文件路径
    %
    % 输出:
    %   无直接输出，但结果将保存到指定的文本文件中
    
    %% 1. 验证输入参数
    % 验证Sigma
    if ~ismatrix(Sigma) || size(Sigma,1) ~= 3 || size(Sigma,2) ~= 3
        error('输入的Sigma必须是一个3x3的矩阵。');
    end
    
    % 验证w
    if ~isvector(w) || length(w) ~= 3
        error('输入的w必须是一个长度为3的向量。');
    end
    w = reshape(w, [3,1]);  % 确保w是3x1列向量
    
    % 验证alpha
    if ~isscalar(alpha) || alpha <= 0 || alpha >= 1
        error('输入的alpha必须是一个介于0和1之间的标量。');
    end
    
    % 验证output_file
    if ~ischar(output_file) && ~isstring(output_file)
        error('输入的output_file必须是一个字符串路径。');
    end
    output_file = char(output_file);  % 确保output_file是字符数组
    
    %% 2. 调整协方差矩阵中的元素
    min_threshold = 1e-9;
    Sigma_adj = Sigma;
    small_elements = Sigma_adj < min_threshold;
    if any(small_elements, 'all')
        Sigma_adj(small_elements) = min_threshold;
        fprintf('已将Sigma中小于 %.1e 的元素替换为 %.1e。\n', min_threshold, min_threshold);
    else
        fprintf('Sigma中所有元素均不小于 %.1e。\n', min_threshold);
    end
    
    %% 3. 计算sigma_L = w' * Sigma * w
    sigma_L = w' * Sigma_adj * w;
    
    %% 4. 计算CVaR
    z_alpha = norminv(alpha);  % 标准正态分布的逆累积分布函数值
    cvar = (normpdf(z_alpha) / (1 - alpha)) * sqrt(sigma_L);
    Cvar = cvar;
    %% 5. 保存结果到文本文件
    try
        fid = fopen(output_file, 'w');
        if fid == -1
            error('无法打开文件 %s 进行写入。', output_file);
        end
        
        fprintf(fid, '条件风险价值（CVaR）计算结果\n');
        fprintf(fid, '===============================\n\n');
        
        fprintf(fid, '协方差矩阵 Sigma (调整后):\n');
        fprintf(fid, '%.12e\t%.12e\t%.12e\n', Sigma_adj.');
        fprintf(fid, '\n');
        
        fprintf(fid, '权重向量 w:\n');
        fprintf(fid, '%.4f\n', w);
        fprintf(fid, '\n');
        
        fprintf(fid, '置信水平 alpha: %.4f\n\n', alpha);
        
        fprintf(fid, '计算结果:\n');
        fprintf(fid, 'sigma_L = w'' * Sigma * w = %.12e\n', sigma_L);
        fprintf(fid, 'CVaR = (pdf(z_alpha) / (1 - alpha)) * sqrt(sigma_L) = %.12e\n', cvar);
        
        fclose(fid);
        fprintf('结果已成功保存到 %s\n', output_file);
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    
    %% 6. 可选：可视化调整后的协方差矩阵
    visualize = false;  % 如果需要可视化，设置为 true
    if visualize
        visualize_covariance(Sigma_adj, '调整后的协方差矩阵 Sigma');
    end
end

%% 辅助函数

function visualize_covariance(Sigma, title_str)
    % visualize_covariance 可视化协方差矩阵
    %
    % 语法:
    %   visualize_covariance(Sigma, title_str)
    %
    % 输入参数:
    %   Sigma     - 3x3 协方差矩阵
    %   title_str - 图像标题字符串
    
    figure;
    imagesc(Sigma);
    colorbar;
    colormap('jet');
    title(title_str, 'FontSize', 14);
    xlabel('变量', 'FontSize', 12);
    ylabel('变量', 'FontSize', 12);
    axis equal tight;
    
    % 在每个元素上显示数值
    textStrings = num2str(Sigma(:), '%.2e');  % 将数值转换为字符串
    textStrings = strtrim(cellstr(textStrings));  % 去除空格
    [x, y] = meshgrid(1:size(Sigma,2), 1:size(Sigma,1));
    hStrings = text(x(:), y(:), textStrings(:), ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'w');
end

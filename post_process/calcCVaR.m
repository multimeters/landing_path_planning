function CVaR = calcCVaR(Sigma, w, alpha)
% calcCVaR 计算在给定协方差矩阵和权重下，正态分布损失函数的CVaR值
%
%   CVaR = calcCVaR(Sigma, w, alpha)
%
% 输入参数：
%   Sigma: 3x3的协方差矩阵
%   w:     3x1的权重向量
%   alpha: 置信水平, 如0.95
%
% 输出参数：
%   CVaR:  对应alpha水平下的条件在险价值
%
% 假设 L = w' * X ，其中 X ~ N(0, Sigma)，则L ~ N(0, w' * Sigma * w)

    % 检查输入参数
    if size(Sigma,1) ~= 3 || size(Sigma,2) ~= 3
        error('Sigma must be a 3x3 matrix.');
    end
    if length(w) ~= 3
        error('w must be a 3x1 vector.');
    end
    if alpha <= 0 || alpha >= 1
        error('alpha must be between 0 and 1.');
    end

    % 计算L的方差
    variance_L = w' * Sigma * w;
    std_L = sqrt(variance_L);

    % 计算VaR时的标准正态分位点 z_alpha
    z_alpha = norminv(alpha);

    % 标准正态密度在z_alpha处的值
    phi_z_alpha = normpdf(z_alpha);

    % 利用解析公式计算CVaR
    CVaR = std_L * (phi_z_alpha / (1 - alpha));
end

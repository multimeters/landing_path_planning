% U船速 10 km/h，即 2.78 m/s
U = 2.78;

% g 重力加速度
g = 9.8;

% h 水深
h = 1.44;

% 遭遇频率系数
k = 1 - U / sqrt(g * h);

% 加载数据
load('../test_data/pixel_1_1506/omega.mat')
load('../test_data/pixel_1000_1000/main_direction.mat')
load('../test_data/pixel_1_1506/energy_density.mat')
%indices = find(energy_density > 0.005);

% omega_e 为遇到频率系数调整后的频率
%omega_e_tmp = k * omega_1024(indices);  % omega_1024 是一个已经加载的变量，包含频率数据
omega_e_tmp = k * omega_1024; 
omega_e = 2 * pi ./ omega_e_tmp;  % 计算周期
% N=512，表示频率谱的长度
N = 512;

% 方向角 μ
angle_degree = max_direction * 180 / pi;  % 转换为度
angle_degree_normalized = mod(angle_degree + 180, 360) - 180;  % 限制在 -180 到 180 之间
mu = abs(angle_degree_normalized);  % 方向角的值
mu = 0 ;
%mu = 89.0;  % 设置一个具体的方向角值
disp(mu);

% 初始化 RAO 的实部和虚部
RAO_3_real = zeros(1, N);
RAO_3_imag = zeros(1, N);
RAO_4_real = zeros(1, N);
RAO_4_imag = zeros(1, N);
RAO_5_real = zeros(1, N);
RAO_5_imag = zeros(1, N);

% 计算每个频率点的 RAO（heave、roll、pitch）
for i = 1:length(omega_e)
    % 查找每个频率点的 RAO
    [RAO_3_real(i), RAO_3_imag(i)] = find_real_im_part(omega_e(i), mu, 'heave');
    [RAO_4_real(i), RAO_4_imag(i)] = find_real_im_part(omega_e(i), mu, 'roll');
    [RAO_5_real(i), RAO_5_imag(i)] = find_real_im_part(omega_e(i), mu, 'pitch');
end

% 将 RAO 转化为复数
RAO_3 = RAO_3_real + 1i * RAO_3_imag;
RAO_4 = RAO_4_real + 1i * RAO_4_imag;
RAO_5 = RAO_5_real + 1i * RAO_5_imag;
% 方差计算公式：Var(x_i) = sum(a_i * (RAO_i(ω) ^ 2))
max3 = 3.0;  % heave 方向上的阈值为3米
max4 = 0.52;  % roll 30度
max5 = 0.52;  % pitch 30度

% 计算每个自由度的方差和协方差
var_x3 = sum(energy_density .* abs(RAO_3).^2) / (max3 * max3); % heave方差
var_x4 = sum(energy_density .* abs(RAO_4).^2) / (max4 * max4); % roll方差
var_x5 = sum(energy_density .* abs(RAO_5).^2) / (max5 * max5); % pitch方差

% 协方差计算公式：Cov(x_i, x_j) = sum(a_i * (RAO_i(ω) * RAO_j*(ω))^*)
cov_x3_x4 = sum(energy_density .* abs((RAO_3 .* conj(RAO_4)))) / (max3 * max4); % heave-roll协方差
cov_x3_x5 = sum(energy_density .* abs((RAO_3 .* conj(RAO_5)))) / (max3 * max5); % heave-pitch协方差
cov_x4_x5 = sum(energy_density .* abs((RAO_4 .* conj(RAO_5)))) / (max4 * max5); % roll-pitch协方差

% 构造协方差矩阵
Sigma = [
    var_x3, 0.00000000001, cov_x3_x5;
    0.00000000001, 0.00000000001, 0.00000000001;
    cov_x3_x5, 0.00000000001, var_x5
];

% 输出协方差矩阵
disp('Covariance Matrix Sigma:');
disp(Sigma);

% 1. 使用协方差矩阵 Sigma
mu = [0, 0, 0];  % 0均值

% 2. 计算 Sigma 的 Cholesky 分解
L = chol(Sigma, 'lower');  % 'lower' 表示返回下三角矩阵

% 3. 生成标准正态分布样本
num_samples = 10000;
z = randn(num_samples, 3);  % 生成 10000 个标准正态分布样本

% 4. 通过 Cholesky 分解转换为具有协方差 Sigma 的样本
samples = z * L + mu;  % 变换样本，使其具有目标协方差

% 计算三个维度的最大值和最小值
max_val = max(abs(samples(:)));  % 取所有样本中绝对值的最大值
min_val = -max_val;  % 对称于原点

% 计算样本的协方差矩阵
cov_matrix = cov(samples);

% 计算协方差矩阵的特征值和特征向量
[EigVec, EigVal] = eig(cov_matrix);

% 特征值的平方根（标准差）
std_devs = sqrt(diag(EigVal));

% 计算95%置信椭球的尺度
confidence_interval = 1.96;  % 对应95%置信度

% 绘制 3D 散点图
figure;
scatter3(samples(:, 1), samples(:, 2), samples(:, 3), 5, 'filled', 'MarkerFaceAlpha', 0.3);

% 保证坐标轴范围相同
xlim([min_val, max_val]);
ylim([min_val, max_val]);
zlim([min_val, max_val]);

% 添加标签和标题
xlabel('Heave');
ylabel('Roll');
zlabel('Pitch');
title('3D Distribution with 95% Confidence Ellipsoid');
grid on;
axis equal;

% 绘制置信椭球
hold on;
% 计算95%置信椭球的参数并绘制
theta = linspace(0, 2*pi, 100);
phi = linspace(0, pi, 100);
[TH, PHI] = meshgrid(theta, phi);

% 将椭球的参数化表示转换为三维坐标
x_ellipsoid = confidence_interval * std_devs(1) * sin(PHI) .* cos(TH);
y_ellipsoid = confidence_interval * std_devs(2) * sin(PHI) .* sin(TH);
z_ellipsoid = confidence_interval * std_devs(3) * cos(PHI);

% 旋转椭球
ellipsoid_points = [x_ellipsoid(:), y_ellipsoid(:), z_ellipsoid(:)] * EigVec';
x_ellipsoid = reshape(ellipsoid_points(:, 1), size(x_ellipsoid));
y_ellipsoid = reshape(ellipsoid_points(:, 2), size(y_ellipsoid));
z_ellipsoid = reshape(ellipsoid_points(:, 3), size(z_ellipsoid));

% 绘制置信椭球
h = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 设置视角和显示比例
view(30, 30);

% 添加图例
hold on;
h_scatter = scatter3(samples(:, 1), samples(:, 2), samples(:, 3), 5, 'filled', 'MarkerFaceAlpha', 0.3);
h_ellipsoid = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 设置图例字体大小
%legend([h_scatter, h_ellipsoid], {'Samples', '95% Confidence Ellipsoid'}, 'Location', 'best', 'FontSize', 14);
legend([h_scatter, h_ellipsoid], {'Samples', '95% Confidence Ellipsoid'}, 'FontSize', 14, 'Position', [0.75, 0.75, 0.2, 0.2]);

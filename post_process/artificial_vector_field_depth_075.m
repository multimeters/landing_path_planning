% MATLAB 程序: visualize_with_resolution.m
% 该程序读取深度数据并绘制带有向量场的图像，正确处理分辨率

% 清理环境
clear; clc; close all;

%% 1. 读取深度数据
depth_file = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/source/depth_a15_512x512.txt';
depth_data = load(depth_file);

% 设置目标深度和阈值
target_depth = 0.0;

% 获取图像的行数和列数
[num_rows, num_cols] = size(depth_data);

% 定义分辨率
dx = 1;  % X方向分辨率（米）
dy = 2;  % Y方向分辨率（米）

% 生成对应的 X 和 Y 网格（以米为单位）
xx = (0:num_cols-1) * dx;  % X轴坐标
yy = (0:num_rows-1) * dy;  % Y轴坐标

[X_depth, Y_depth] = meshgrid(xx, yy);  % 创建网格

%% 2. 初始化并生成“人工”向量场（基于深度数据）
U = zeros(num_rows, num_cols);  % 水平方向上的向量分量
V = zeros(num_rows, num_cols);  % 垂直方向上的向量分量

% 设置最大值和最小值
max_length = 2;
min_length = 0;

% 缩放因子，用来调整向量的长度
scale_factor = 0.8;  % 可以根据需要调整

% 对每一行，找到最接近目标深度的位置，并构建向量场
for i = 1:num_rows
    for j = 1:num_cols
        % 如果深度大于 0.0
        if depth_data(i, j) > target_depth
            % 计算向量长度：深度越小，向量长度越长
            vector_length = scale_factor / (depth_data(i, j) - target_depth);
            
            % 限制向量长度在 [min_length, max_length] 范围内
            vector_length = max(min(vector_length, max_length), min_length);
            
            % 向量方向是水平向左，所以水平方向分量 U 为负数
            U(i, j) = -vector_length;
            V(i, j) =  0;
        end
    end
end

%% 3. 读取已存在的向量场（原始向量场）
u_vf_file = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/u_vf.txt';
v_vf_file = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/v_vf.txt';
U_vf = load(u_vf_file);
V_vf = load(v_vf_file);

% 确保两个向量场大小一致
if any(size(U_vf) ~= size(U))  % 用 any(...) 避免逻辑误判
    error('生成的向量场和已加载的向量场大小不一致！');
end

% 叠加向量场
U_combined = U + U_vf;
V_combined = V + V_vf;

%% 4. 可视化
% 4.1 显示深度图和“原始向量场”
figure;  % 打开新窗口
imagesc(xx, yy, depth_data);  % 使用 imagesc 显示深度图，坐标为米单位
hold on;

% 添加热力图的颜色条
colorbar;
colormap jet;  % 使用 Jet 色图表示深度
caxis([min(depth_data(:)), max(depth_data(:))]);  % 设置颜色条范围

% 设置坐标轴范围（X轴和Y轴使用米为单位）
axis([xx(1) xx(end) yy(1) yy(end)]);  
axis equal;  
set(gca, 'YDir', 'normal'); 
xlabel('X (米)');
ylabel('Y (米)');
title('深度热力图 + 原始向量场');

% 绘制“原始”向量场（U_vf, V_vf）
% 注意：V_vf * dy，目的是在可视化时根据实际 dy 拉伸/缩放垂直方向的箭头
quiver(X_depth, Y_depth, U_vf, V_vf * dy, 'Color', 'w', 'MaxHeadSize', 2);
hold off;

% 4.2 显示深度图和“人工向量场”（即刚才在第2步中根据深度生成的 U 和 V）
figure;  % 打开新窗口
imagesc(xx, yy, depth_data);
hold on;

colorbar;
colormap jet;
caxis([min(depth_data(:)), max(depth_data(:))]);

axis([xx(1) xx(end) yy(1) yy(end)]);
axis equal;
set(gca, 'YDir', 'normal');
xlabel('X (米)');
ylabel('Y (米)');
title('深度热力图 + 人工向量场');

% 绘制人工向量场（U, V）
quiver(X_depth, Y_depth, U, V * dy, 'Color', 'w', 'MaxHeadSize', 2);
hold off;

% 4.3 显示深度图和“合成向量场”（U_combined, V_combined）
figure;  % 打开新窗口
imagesc(xx, yy, depth_data);
hold on;

colorbar;
colormap jet;
caxis([min(depth_data(:)), max(depth_data(:))]);

axis([xx(1) xx(end) yy(1) yy(end)]);
axis equal;
set(gca, 'YDir', 'normal');
xlabel('X (米)');
ylabel('Y (米)');
title('深度热力图 + 合成向量场');

% 绘制合成向量场
quiver(X_depth, Y_depth, U_combined, V_combined * dy, 'Color', 'w', 'MaxHeadSize', 2);
hold off;

%% 5. 保存合成向量场为文本文件
u_combined_file = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/u_combined.txt';
v_combined_file = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/vf/v_combined.txt';

save(u_combined_file, 'U_combined', '-ascii');
save(v_combined_file, 'V_combined', '-ascii');

disp('合成向量场已保存！');

% MATLAB 程序：绘制 direction_differences_with_cvar.txt 文件的第6列数据

% 清空工作区和命令窗口
clear;
clc;

% 指定文件名（确保文件位于当前工作目录或提供完整路径）
filename = 'output_data_path_e_01_lam_1001_sam_1000000/direction_differences_with_cvar_path_37.txt';

% 检查文件是否存在
if ~isfile(filename)
    error('文件 %s 不存在。请检查文件名和路径。', filename);
end

% 加载数据
% 假设数据以空格、制表符或逗号分隔
data = load(filename);

% 检查数据是否至少有6列
[numRows, numCols] = size(data);
if numCols < 6
    error('数据文件中列数不足6列。当前列数：%d。', numCols);
end

% 提取第6列数据
column6 = data(:,6);

% 绘制第6列数据
figure;
plot(column6, 'b-', 'LineWidth', 1.5); % 使用蓝色实线，线宽1.5
xlabel('样本索引', 'FontSize', 12);
ylabel('第6列数值', 'FontSize', 12);
title('direction\_differences\_with\_cvar.txt 第6列数据', 'FontSize', 14);
grid on; % 添加网格
legend('第6列数据', 'Location', 'best');

% 可选：保存图形为PNG文件
% saveas(gcf, 'column6_plot.png');

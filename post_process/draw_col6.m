% 1. 清空工作区
clear; clc; close all;

% 2. 加载数据
data = load('\\wsl.localhost\ubuntu\home\lhl\share\rip\scripts\landing_path_planning\post_process\output_data_path_e_01_lam_1001_sam_1000000\direction_differences_with_cvar.txt');

% 3. 提取第 6 列数据（MATLAB 中索引从 1 开始，第 6 列在代码中要用 data(:,6)）
% 如果你的文件实际需要的是第 5 列索引位置，请根据需求改为 data(:,5)
col6 = data(:,5);  

% 4. 绘制曲线
figure;
p = plot(col6, '-o', 'LineWidth', 1.5); 
xlabel('索引');  
ylabel('数据值');
title('第 6 列数据曲线');
grid on;

% 5. 打开数据光标模式（Data Cursor），并设置自定义回调函数
dcmObj = datacursormode(gcf);
set(dcmObj, 'UpdateFcn', @(src, event) myCallbackFcn(src, event, data));

%% ---------------------------- 回调函数 --------------------------------
function output_txt = myCallbackFcn(~, event, data)
    % 当鼠标点击曲线某一点时，返回该点对应的信息
    
    % 获取点击的点在数据中的索引
    idx = event.DataIndex;  
    
    % 这里演示输出 data(:,1) 和 data(:,2)，当然也可以加更多需要显示的内容
    x_val = data(idx, 1);
    y_val = data(idx, 2);
    
    % 自定义要显示在光标提示框中的文本
    output_txt = { ...
        ['Index: ', num2str(idx)], ...
        ['data(:,1) = ', num2str(x_val)], ...
        ['data(:,2) = ', num2str(y_val)] ...
        };
end

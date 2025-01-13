% batch_process_paths.m
% 
% 功能:
% 该脚本遍历指定文件夹中的所有路径文件，并对每个文件调用
% combined_plot_and_process_path 函数进行处理。
%
% 所有输出文件（包括 direction_differences.txt）将保存在统一的
% 'output_data' 文件夹中，文件名包含路径文件的基名以确保唯一性。
%
% 使用方法:
% 在 MATLAB 命令窗口中运行此脚本，确保 combined_plot_and_process_path.m
% 文件在当前工作目录或 MATLAB 路径中。

%% ---------------------------- 初始化部分 ----------------------------

% 清空工作区和图形窗口
clear;
clc;
close all;

%% --------------------- 设置参数 ---------------------

% 定义包含所有路径文件的文件夹路径
paths_folder = '\\wsl.localhost\Ubuntu\home\lhl\share\rip\scripts\landing_path_planning\VFRRT_Planner\path_e_01_lam_1001_sam_1000000\';

% 定义路径文件的命名模式，例如 'path_*.txt'
path_file_pattern = 'path_*.txt';

% 定义 main_directions.mat 文件的路径
main_directions_path = '/home/lhl/share/rip/scripts/landing_path_planning/vf_estimate/data/output/main_directions.mat';

% 定义主输出文件夹（统一输出到 paths_folder 下的 'output_data' 文件夹）
master_output_folder = 'output_data_path_e_01_lam_1001_sam_1000000';

% 定义采样间隔（米）
sample_interval = 1;

% 定义绘图保存的文件名模板，例如 'path_with_tangents_<path_name>.png'
figure_filename_template = 'path_with_tangents_%s.png'; % %s 将被替换为路径文件的基名

% 定义 x 和 y 坐标转换为索引的缩放因子
scale_x = 1;
scale_y = 2;

%% --------------------- 创建主输出文件夹 ---------------------

% 创建主输出文件夹（如果不存在）
if ~exist(master_output_folder, 'dir')
    mkdir(master_output_folder);
    fprintf('已创建主输出文件夹: %s。\n', master_output_folder);
else
    fprintf('主输出文件夹已存在: %s。\n', master_output_folder);
end

%% --------------------- 获取所有路径文件列表 ---------------------

% 构建完整的文件搜索模式
search_pattern = fullfile(paths_folder, path_file_pattern);

% 获取所有匹配的路径文件
path_files = dir(search_pattern);

% 检查是否找到任何文件
if isempty(path_files)
    error('在文件夹 %s 中未找到匹配的路径文件。', paths_folder);
end

fprintf('在文件夹 %s 中找到 %d 个路径文件。\n', paths_folder, length(path_files));

%% --------------------- 批量处理每个路径文件 ---------------------

for k = 1:length(path_files)
    % 获取当前路径文件的名称和完整路径
    current_path_file = path_files(k).name;
    current_path_full = fullfile(paths_folder, current_path_file);
    
    % 提取路径文件的基名（不带扩展名）
    [~, path_name, ~] = fileparts(current_path_file);
    
    % 定义采样点输出文件的路径
    % 将采样点文件保存在主输出文件夹下的 'sampled_<path_name>.txt'
    sampled_path_filename = sprintf('sampled_path_%s.txt', path_name);
    sampled_path_full = fullfile(master_output_folder, sampled_path_filename);
    
    % 定义绘图文件名，例如 'path_with_tangents_<path_name>.png'
    figure_filename = sprintf(figure_filename_template, path_name);
    
    % 定义切向量夹角文件名，例如 'tangent_angles_<path_name>.txt'
    angles_filename = sprintf('tangent_angles_%s.txt', path_name);
    
    % 定义方向差异文件名，例如 'direction_differences_<path_name>.txt'
    direction_differences_filename = sprintf('direction_differences_%s.txt', path_name);
    
    %% --------------------- 调用处理函数 ---------------------
    
    fprintf('\n===== 处理文件 %s (%d/%d) =====\n', current_path_file, k, length(path_files));
    try
        combined_plot_and_process_path(...
            current_path_full, ...                    % original_path_file
            sampled_path_full, ...                    % sampled_path_file
            main_directions_path, ...                 % main_directions_path
            'OutputFolder', master_output_folder, ... % 输出文件夹为统一的主输出文件夹
            'SampleInterval', sample_interval, ...    % 采样间隔
            'FigureFilename', figure_filename, ...    % 绘图文件名
            'AnglesFilename', angles_filename, ...    % 切向量夹角文件名
            'DirectionDifferencesFilename', direction_differences_filename, ... % 方向差异文件名
            'ScaleX', scale_x, ...                    % x 轴缩放因子
            'ScaleY', scale_y ...                     % y 轴缩放因子
        );
        fprintf('成功处理文件: %s。\n', current_path_file);
    catch ME
        warning('处理文件 %s 时出错: %s', current_path_file, ME.message);
    end
end

fprintf('\n===== 所有路径文件已处理完毕 =====\n');

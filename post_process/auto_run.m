% 设置不同的 exploration 参数值
exploration_values = 0.1;
initial_lambda = [2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0];
initial_lambda_samples = 10000;

% MATLAB 脚本路径
matlab_script = 'plot_v1.m';  % 这个脚本不再使用了

% 循环执行命令
for i = 1:length(initial_lambda)
    initial_lambda_i = initial_lambda(i)+16.0;
    
    % 构建命令行来执行 ./your_program
    command = sprintf('wsl.exe /home/lhl/share/rip/scripts/landing_path_planning/VFRRT_Planner/build/VFRRTPlanner -e %.1f -l %.1f -s %.0f', ...
                      exploration_values, initial_lambda_i, initial_lambda_samples);

    % 执行 ./your_program
    fprintf('执行: %s\n', command);
    [status, result] = system(command);
    
    % 检查命令执行的状态
    if status == 0
        fprintf('程序成功执行，输出：\n%s\n', result);
    else
        fprintf('程序执行失败，错误信息：\n%s\n', result);
    end
    
    % 构建保存结果的文件名
    outputPath = sprintf('path_planning_exploration_%.1f_lambda_%.1f_samples_%.0f.png', exploration_values, initial_lambda_i, initial_lambda_samples);
    
    % 调用 visualize_all_paths_with_vector_field_and_depth(outputPath) 并传递文件名
    fprintf('执行 visualize_all_paths_with_vector_field_and_depth 函数, 保存路径: %s\n', outputPath);
    visualize_all_paths_with_vector_field_and_depth(outputPath);
    
    % 可选：根据需要添加延时
    pause(1);  % 延时 1 秒，可以根据需求调整
end

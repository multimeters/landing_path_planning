% 设置文件路径
original_dir = ''; % 修改为原始txt文件所在的目录
mat_dir_eta = 'pixel_data/eta';      % 修改为eta.mat文件所在的目录
mat_dir_u = 'pixel_data/u';          % 修改为u.mat文件所在的目录
mat_dir_v = 'pixel_data/v';          % 修改为v.mat文件所在的目录

% 定义文件数量和数据尺寸
num_files = 1000;
rows = 4096;
cols = 2048;

% 采样次数
num_samples = 100;

% 初始化随机数种子以确保可重复性
rng(42);

% 随机生成采样索引
sample_files = randi([0 num_files-1], [num_samples, 1]);
sample_rows = randi([1 rows], [num_samples, 1]);
sample_cols = randi([1 cols], [num_samples, 1]);

% 初始化错误计数
error_count_eta = 0;
error_count_u = 0;
error_count_v = 0;


% 对eta进行测试
error_count_eta = test_variable('eta', original_dir, mat_dir_eta, sample_files, sample_rows, sample_cols, num_samples);

% 对u进行测试
error_count_u = test_variable('u', original_dir, mat_dir_u, sample_files, sample_rows, sample_cols, num_samples);

% 对v进行测试
error_count_v = test_variable('v', original_dir, mat_dir_v, sample_files, sample_rows, sample_cols, num_samples);

% 输出测试结果
if error_count_eta == 0
    fprintf('All sampled values match between .txt and .mat files for eta.\n');
else
    fprintf('Total mismatches for eta: %d\n', error_count_eta);
end

if error_count_u == 0
    fprintf('All sampled values match between .txt and .mat files for u.\n');
else
    fprintf('Total mismatches for u: %d\n', error_count_u);
end

if error_count_v == 0
    fprintf('All sampled values match between .txt and .mat files for v.\n');
else
    fprintf('Total mismatches for v: %d\n', error_count_v);
end

% 定义函数来测试不同的变量
function error_count = test_variable(var_name, original_dir, mat_dir, sample_files, sample_rows, sample_cols, num_samples)
    error_count = 0;
    for i = 1:num_samples
        file_index = sample_files(i);
        row_index = sample_rows(i);
        col_index = sample_cols(i);

        % 构造文件名
        original_file = fullfile(original_dir, sprintf('%s_%05d', var_name, file_index));
        mat_file = fullfile(mat_dir, sprintf('%s_%05d.mat', var_name, file_index));

        % 从txt文件中读取相应的值
        fid = fopen(original_file, 'r');
        data_txt = fscanf(fid, '%f', [cols, rows])';
        fclose(fid);

        value_txt = data_txt(row_index, col_index);

        % 从.mat文件中读取相应的值
        data_mat = load(mat_file);
        value_mat = data_mat.data(row_index, col_index);

        % 检查是否一致
        if value_txt ~= value_mat
            fprintf('Mismatch in file %s_%05d.mat at row %d, col %d: txt value = %f, mat value = %f\n', ...
                var_name, file_index, row_index, col_index, value_txt, value_mat);
            error_count = error_count + 1;
        end
    end
end

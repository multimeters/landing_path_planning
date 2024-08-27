function convert_txt_to_mat(prefix, index)
    % 生成文件名
    filename = sprintf('%s_%05d', prefix, index);
    
    % 读取数据，指定空格为分隔符
    data = importdata(filename, ' ');
    
    % 如果数据是结构体，从结构体的数据字段中提取数值
    if isstruct(data)
        data = data.data;
    end
    
    % 创建保存路径
    save_path = fullfile('pixel_data', prefix);
    
    % 如果文件夹不存在，则创建
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    % 保存为 .mat 文件
    save(fullfile(save_path, sprintf('%s_%05d.mat', prefix, index)), 'data');
end

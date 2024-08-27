delete(gcp('nocreate')); % 删除现有的并行池（如果存在）
parpool; % 启动并行池

parfor k = 0:999
    k
    convert_txt_to_mat('eta', k);
end

parfor i = 0:999
    i
    convert_txt_to_mat('u', i);
end

parfor j = 0:999
    j
    convert_txt_to_mat('v', j);
end

delete(gcp('nocreate')); % 关闭并行池

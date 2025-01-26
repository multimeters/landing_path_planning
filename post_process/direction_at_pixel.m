% 假设 energy_distribution 是一个 72x1 的变量
% 角度范围从 -180 到 180
angles = linspace(-180, 180, 72); % 创建一个角度范围从 -180 到 180 的向量

% 将角度转换为弧度
angles_rad = deg2rad(angles);

% 计算加权平均角度
% 使用能量作为权重，计算每个角度的加权平均
weighted_angle = sum(energy_distribution .* cos(angles_rad)) / sum(energy_distribution);
weighted_angle_sin = sum(energy_distribution .* sin(angles_rad)) / sum(energy_distribution);

% 计算主方向的角度（反正切）
main_direction_rad = atan2(weighted_angle_sin, weighted_angle);

% 将主方向从弧度转换回度
main_direction_deg = rad2deg(main_direction_rad);

% 确保主方向在 -180 到 180 度之间
if main_direction_deg > 180
    main_direction_deg = main_direction_deg - 360;
elseif main_direction_deg < -180
    main_direction_deg = main_direction_deg + 360;
end

% 显示主方向
disp(['主方向: ', num2str(main_direction_deg), '°']);

% 绘制能量分布图
figure;
polarplot(angles_rad, energy_distribution, '-o'); % 极坐标图
hold on;

% 标示主方向
polarplot(main_direction_rad, max(energy_distribution), 'r*', 'MarkerSize', 10);

% 标题
title(['主方向: ', num2str(main_direction_deg), '°']);

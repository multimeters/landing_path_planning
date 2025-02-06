% 假设数据已经存储在一个名为 data 的向量中
data = randn(1, 1024);  % 示例数据，可以根据需要替换成你的实际数据

% 设置采样频率
Fs = 1;  % 采样频率为1 Hz

% 计算FFT
N = length(data);  % 数据长度（1024）
Y = fft(data, N);   % 计算FFT

% 计算频率轴
f = (0:N/2-1)*(Fs/N);  % 频率分量，从0到采样频率的一半
omega = 2*pi*f;        % 转换为角频率（弧度/秒）
Y_half = Y(1:N/2);     % FFT结果的前半部分（对称）

% 计算前512个频率和对应的FFT幅度
omega_1024 = omega(1:512);  % 角频率的前512个
Y_512 = abs(Y_half(1:512));  % 对应的FFT幅度（取前512）

% 绘制频谱
figure;
plot(omega_1024, Y_512);
title('前512个频率的FFT幅度 (采样频率 1Hz, 角频率 rad/s)');
xlabel('角频率 (\omega) [rad/s]');
ylabel('幅度');
grid on;

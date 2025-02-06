import numpy as np
import time
import os
from scipy.io import loadmat, savemat

# 基本参数
fs = 1
N_1024 = 1024
rowsPerSegment = 256
numSegments = 1

# 频率与角频率向量
f_1024 = fs * np.arange(N_1024 / 2) / N_1024
f_1024 = f_1024.astype(float)
omega_1024 = 2 * np.pi * f_1024
if len(omega_1024) > 1:
    delta_omega = omega_1024[1] - omega_1024[0]
    print(delta_omega)
else:
    delta_omega = 0

output_dir = '../data/processing/fft_results'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for i in range(numSegments):
    rowStart = i * rowsPerSegment
    rowEnd = (i + 1) * rowsPerSegment

    eta_filename = f"../data/processing/e_output/eta_0-255.mat"
    mask_filename = f"../data/processing/mask/mask.mat"

    # 读取MAT文件（CPU内存）
    eta_data = loadmat(eta_filename)
    eta_block_cpu = eta_data['data']  # shape: (256, 2048, 1024)

    mask_data = loadmat(mask_filename)
    mask_block_cpu = mask_data['mask']  # shape: (256, 2048)

    # 利用mask在CPU上处理：对mask为0的位置，将eta数据全置0
    mask_bool = (mask_block_cpu == 0)  # shape: (256, 2048)
    # 扩展mask以匹配eta_block的形状
    mask_bool_expanded = mask_bool[..., np.newaxis]  # shape: (256, 2048, 1)
    # 利用广播将mask应用到eta_block
    #eta_block_cpu = eta_block_cpu * mask_bool_expanded  # 直接将不需要处理的部分置0

    # 开始FFT
    start_time = time.time()
    fft_wave_height_block = np.fft.fft(eta_block_cpu, axis=2)
    fft_time_1024 = time.time() - start_time
    print(f"第{i + 1}块FFT计算时间（1024个数据，CPU版）：{fft_time_1024:.6f}秒")

    # 计算能量密度
    abs_fft = np.abs(fft_wave_height_block)
    energy_density_block_tmp = (abs_fft ** 2) / N_1024
    energy_density_block = energy_density_block_tmp * delta_omega

    # 只取正半部分频率
    fft_half = fft_wave_height_block[:, :, :N_1024 // 2]
    energy_density_half = energy_density_block[:, :, :N_1024 // 2]

    # 将结果保存为MAT文件
    fft_output_filename = os.path.join(output_dir, f"fft_block_{rowStart}-{rowEnd - 1}.mat")
    energy_output_filename = os.path.join(output_dir, f"energy_block_{rowStart}-{rowEnd - 1}.mat")

    savemat(fft_output_filename, {'fft_half': fft_half})
    savemat(energy_output_filename, {'energy_density_half': energy_density_half})

print("所有块的FFT处理完成。")

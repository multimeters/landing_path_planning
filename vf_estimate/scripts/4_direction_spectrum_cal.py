import numpy as np
from scipy.io import loadmat, savemat
import glob
import time
import sys
from numba import njit, prange

# ----------------------------- 设置 Matplotlib 后端 -----------------------------
import matplotlib
matplotlib.use('Agg')  # 使用 'Agg' 后端以避免 Qt 相关错误
import matplotlib.pyplot as plt

# ----------------------------- 配置 -----------------------------

# 定义方向分箱的数量
NUM_BINS = 36

# 定义用于绘图的特定文件索引和像素坐标
# 根据您的数据和需求调整这些值
TARGET_FILE_INDEX = 0  # 从0开始的索引（例如，4对应第5个文件）
TARGET_ROW = 206        # 从0开始的索引
TARGET_COL = 59         # 从0开始的索引

# 定义块大小（每个块的行数）
# 根据系统的内存容量进行调整
CHUNK_SIZE = 64

# ----------------------------- 开始计时 -----------------------------

start_time = time.time()

# ----------------------------- 文件处理 -----------------------------

# 检索并排序所有相关的 .mat 文件
data_files = sorted(glob.glob('../data/processing/e_output/eta_*.mat'))
angle_files = sorted(glob.glob('../data/processing/a_output/a_*.mat'))

# 验证是否找到必要的文件
if not data_files:
    print("错误: 未找到任何 'e_output/eta_*.mat' 文件。请检查文件路径和文件是否存在。")
    sys.exit(1)

if not angle_files:
    print("错误: 未找到任何 'a_output/a_*.mat' 文件。请检查文件路径和文件是否存在。")
    sys.exit(1)

# 确保数据文件的数量与角度文件的数量匹配
if len(data_files) != len(angle_files):
    print(f"错误: 'e_output/eta_*.mat' 文件数量 ({len(data_files)}) 与 'a_output/a_*.mat' 文件数量 ({len(angle_files)}) 不匹配。")
    sys.exit(1)

# ----------------------------- 定义方向分箱 -----------------------------

# 从 -π 到 π 创建线性间隔的方向分箱
direction_bins = np.linspace(-np.pi, np.pi, NUM_BINS + 1, dtype=np.float32)  # 37 个边缘
bin_centers = (direction_bins[:-1] + direction_bins[1:]) / 2  # 36 个中心

# 确保数据类型与 Numba 兼容
direction_bins = direction_bins.astype(np.float32)
bin_centers = bin_centers.astype(np.float32)

# ----------------------------- Numba 加速处理 -----------------------------

@njit(parallel=True)
def process_file_numba(e_data, a_data, direction_bins, bin_centers, num_bins):
    """
    处理单个文件的数据以计算主要方向和能量分布。

    参数:
    - e_data: (rows, cols, depth) float32 数组，表示波高。
    - a_data: (rows, cols, depth) float32 数组，表示角度。
    - direction_bins: (num_bins +1,) float32 数组，定义方向分箱的边缘。
    - bin_centers: (num_bins,) float32 数组，定义方向分箱的中心。
    - num_bins: 整数，表示方向分箱的数量。

    返回:
    - main_dirs: (rows, cols) float32 数组，每个像素的主要方向。
    - energy_sum: (rows, cols, num_bins) float32 数组，归一化的能量分布。
    """
    rows, cols, depth = e_data.shape
    main_dirs = np.zeros((rows, cols), dtype=np.float32)
    energy_sum = np.zeros((rows, cols, num_bins), dtype=np.float32)
    
    for i in prange(rows):
        for j in prange(cols):
            # 找到当前像素的最小波高
            min_height = e_data[i, j, 0]
            for k in range(1, depth):
                if e_data[i, j, k] < min_height:
                    min_height = e_data[i, j, k]
            
            # 调整波高并计算能量
            total_energy = 0.0
            energies = np.zeros(depth, dtype=np.float32)
            for k in range(depth):
                adjusted = e_data[i, j, k] + np.abs(min_height)
                energies[k] = adjusted * adjusted
                total_energy += energies[k]
            
            # 归一化能量
            for k in range(depth):
                energies[k] /= total_energy
            
            # 将能量分配到相应的方向分箱
            for k in range(depth):
                angle = a_data[i, j, k]
                bin_idx = -1
                # 手动分箱，因为 np.digitize 不被 Numba 支持
                for b in range(num_bins):
                    if direction_bins[b] <= angle < direction_bins[b + 1]:
                        bin_idx = b
                        break
                if bin_idx == -1:
                    if angle == direction_bins[-1]:
                        bin_idx = num_bins -1
                    else:
                        bin_idx = 0  # 超出范围时的回退
                # 累积能量
                energy_sum[i, j, bin_idx] += energies[k]
            
            # 基于最大能量确定主要方向
            max_bin = 0
            max_energy = energy_sum[i, j, 0]
            for b in range(1, num_bins):
                if energy_sum[i, j, b] > max_energy:
                    max_energy = energy_sum[i, j, b]
                    max_bin = b
            main_dirs[i, j] = bin_centers[max_bin]
    
    return main_dirs, energy_sum

# ----------------------------- 初始化结果存储 -----------------------------

# 根据第一个文件确定数据的形状
try:
    first_e_data = loadmat(data_files[0])['data']
    rows_per_file, cols_per_file, depth = first_e_data.shape
except KeyError:
    print(f"错误: 文件 {data_files[0]} 中未找到关键字 'data'。")
    sys.exit(1)
except IndexError:
    print(f"错误: 文件 {data_files[0]} 中的 'data' 形状不正确。")
    sys.exit(1)

# 根据文件数量和每个文件的行数计算总行数
total_files = len(data_files)
total_rows = rows_per_file * total_files
total_cols = cols_per_file  # 假设所有文件具有相同的列数

# 初始化 main_directions 数组
# 使用 float32 来平衡精度和内存使用
main_directions = np.zeros((total_rows, total_cols), dtype=np.float32)

# ----------------------------- 处理循环 -----------------------------

# 初始化用于跟踪进度和绘图的变量
image_index = 1
current_row = 0

for file_index in range(len(data_files)):
    data_file = data_files[file_index]
    angle_file = angle_files[file_index]
    print(f"Processing file {file_index+1}/{len(data_files)}: {data_file}")

    # ----------------------------- 加载数据 -----------------------------

    # 加载波高数据 (e_data)
    try:
        e_data = loadmat(data_file)['data'].astype(np.float32)  # 形状: (rows, cols, depth)
    except KeyError:
        print(f"错误: 文件 {data_file} 中未找到关键字 'data'。")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 无法读取文件 {data_file}。详细信息: {e}")
        sys.exit(1)

    # 加载角度数据 (a_data)
    try:
        a_data = loadmat(angle_file)['data'].astype(np.float32)  # 形状: (rows, cols, depth)
    except KeyError:
        print(f"错误: 文件 {angle_file} 中未找到关键字 'data'。")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 无法读取文件 {angle_file}。详细信息: {e}")
        sys.exit(1)

    # ----------------------------- 处理数据 -----------------------------

    # 使用 Numba 加速的函数处理当前文件
    main_dirs, energy_distribution = process_file_numba(e_data, a_data, direction_bins, bin_centers, NUM_BINS)

    # 将处理后的主要方向分配到 main_directions 数组
    main_directions[current_row:current_row + main_dirs.shape[0], :] = main_dirs

    # ----------------------------- 绘图 -----------------------------

    # 为目标文件中的特定像素生成绘图
    # 根据需要调整 TARGET_FILE_INDEX, TARGET_ROW, 和 TARGET_COL
    print(file_index)
    print(TARGET_FILE_INDEX)
    if file_index == TARGET_FILE_INDEX:
        # 确保目标行和列在范围内
        if TARGET_ROW >= energy_distribution.shape[0] or TARGET_COL >= energy_distribution.shape[1]:
            print(f"警告: 目标行或列超出范围。跳过绘图。")
        else:
            energy_dist = energy_distribution[TARGET_ROW, TARGET_COL, :]
            main_dir = main_dirs[TARGET_ROW, TARGET_COL]

            plt.figure(figsize=(10, 6))
            plt.bar(np.degrees(bin_centers), energy_dist, width=10, color='skyblue', edgecolor='black')
            plt.title(f'Energy Distribution and Main Direction (Row {TARGET_ROW+1}, Column {TARGET_COL+1})')
            plt.xlabel('Direction (degrees)')
            plt.ylabel('Normalized Energy')
            plt.axvline(np.degrees(main_dir), color='red', linestyle='--', label=f'Main Direction: {np.degrees(main_dir):.2f}°')
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            plt.savefig(f'../data/output/distribution_{image_index}_row_{TARGET_ROW+1}_col_{TARGET_COL+1}.png')
            plt.close()

            print(f"已保存能量分布图像: distribution_{image_index}_row_{TARGET_ROW+1}_col_{TARGET_COL+1}.png")
            image_index += 1

    # ----------------------------- 清理 -----------------------------

    # 更新 main_directions 数组的当前行指针
    current_row += main_dirs.shape[0]

    # 通过删除大型变量释放内存
    del e_data, a_data, main_dirs, energy_distribution

    # 可选: 强制进行垃圾回收（如有必要，请取消注释）
    # import gc
    # gc.collect()

# ----------------------------- 保存结果 -----------------------------

# 将 main_directions 数组保存到 .mat 文件中
savemat('../data/output/main_directions.mat', {'main_directions': main_directions})
print("已成功保存主方向数据到 'main_directions.mat'。")

# ----------------------------- 计算并保存 u_main 和 v_main -----------------------------

# 计算水平（u）和垂直（v）分量
u_main = np.cos(main_directions)
v_main = np.sin(main_directions)

# 保存 u_main 到 u_main.txt
# 使用逗号分隔，并保留6位小数
np.savetxt('u_main.txt', u_main, fmt='%.6f', delimiter=',')
print("已成功保存水平分量数据到 'u_main.txt'。")

# 保存 v_main 到 v_main.txt
# 使用逗号分隔，并保留6位小数
np.savetxt('v_main.txt', v_main, fmt='%.6f', delimiter=',')
print("已成功保存垂直分量数据到 'v_main.txt'。")

# ----------------------------- 结束计时 -----------------------------

end_time = time.time()
total_time = end_time - start_time
print(f"总执行时间: {total_time:.2f} 秒")
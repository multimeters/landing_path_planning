import numpy as np
import scipy.io as sio
import os
from concurrent.futures import ProcessPoolExecutor

def process_single_file(i, per_file_template, rao_data, max_row, step, per_dir, output_dir):
    """
    处理单个 per 文件并保存结果
    """
    start = i * step
    end = (i + 1) * step - 1
    per_mat_path = os.path.join(per_dir, f'per_{start}-{end}.mat')  # 从 per_dir 文件夹中读取文件
    output_mat_path = os.path.join(output_dir, f'heave_rao-{start}-{end}.mat')
    
    print(f"\n处理文件: {per_mat_path} ...")
    
    # 检查 per 文件是否存在
    if not os.path.isfile(per_mat_path):
        print(f"文件 {per_mat_path} 不存在，跳过...")
        return
    
    # 加载 per 文件
    per_mat = sio.loadmat(per_mat_path)
    
    # 假设变量名为 'per'，请根据实际情况修改
    if 'per' in per_mat:
        per_data = per_mat['per']
    else:
        print(f"{per_mat_path} 中的变量名有:", per_mat.keys())
        raise KeyError(f"未找到变量名 'per' 在 {per_mat_path} 文件中。请检查变量名。")
    
    # 确认 per_data 的维度为 (256, 512, 512)
    expected_shape = (step, 512, 512)
    if per_data.shape != expected_shape:
        raise ValueError(f"per_data 的形状不符合预期 {expected_shape}，实际形状为 {per_data.shape}")
    
    # 计算 row_index
    print("计算 row_index...")
    row_index = np.floor(per_data / 0.1).astype(int)
    row_index = row_index - 1
    # 确保 row_index 在有效范围内
    print("确保 row_index 在有效范围内...")
    row_index = np.clip(row_index, 0, max_row)
    
    # 提取 heave 数据
    print("提取 heave 数据...")
    heave_data = rao_data[row_index, 3]
    
    # 保存为新的 .mat 文件，变量名为 'heave_rao'
    print(f"保存结果到 {output_mat_path} ...")
    sio.savemat(output_mat_path, {'heave_rao': heave_data})
    
    print(f"文件 {output_mat_path} 处理完成！")

def process_multiple_mat_files(per_file_template, rao_mat_path, output_dir, per_dir, num_files=16, step=256):
    """
    读取多个 per_*.mat 文件和 rao_heave_at_0.mat 文件，处理数据并保存为 heave_rao-*.mat 文件。
    
    参数:
    - per_file_template: str, per 文件名模板，例如 'per_{}-{}.mat'
    - rao_mat_path: str, rao_heave_at_0.mat 文件路径
    - output_dir: str, 输出文件的目录
    - per_dir: str, 存放 per 文件的目录
    - num_files: int, 需要处理的文件数量（默认 16）
    - step: int, 每个文件覆盖的范围步长（默认 256）
    """
    
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 加载 rao_heave_at_0.mat 文件一次，以提高效率
    print("加载 rao_heave_at_0.mat 文件...")
    rao_mat = sio.loadmat(rao_mat_path)
    
    # 假设变量名为 'heave'，请根据实际情况修改
    if 'heave' in rao_mat:
        rao_data = rao_mat['heave']
    else:
        print("rao_heave_at_0.mat 中的变量名有:", rao_mat.keys())
        raise KeyError("未找到变量名 'heave' 在 rao_heave_at_0.mat 文件中。请检查变量名。")
    
    print(f"rao_data 的形状: {rao_data.shape}")
    
    # 确保 rao_data 至少有 4 列
    if rao_data.shape[1] < 4:
        raise ValueError(f"rao_data 的列数不足 4，实际列数为 {rao_data.shape[1]}")
    
    # 获取 rao_data 的行数
    max_row = rao_data.shape[0] - 1
    
    # 使用并行处理来加速文件处理
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_single_file, i, per_file_template, rao_data, max_row, step, per_dir, output_dir) 
                   for i in range(num_files)]
        
        # 等待所有任务完成
        for future in futures:
            future.result()
    
    print("\n所有文件处理完成！")

if __name__ == "__main__":
    # 定义参数
    per_file_template = 'per_{}-{}.mat'  # 示例模板，可根据需要调整
    rao_mat_path = '../data/source/rao_heave_at_0.mat'
    output_dir = '../data/processing/output_heave_rao'  # 输出文件夹
    per_dir = '../data/processing/per'  # 存放 per 文件的目录
    num_files = 1
    step = 256
    
    # 调用处理函数
    process_multiple_mat_files(per_file_template, rao_mat_path, output_dir, per_dir, num_files, step)

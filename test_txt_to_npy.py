import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor

# 设置文件路径
original_dir = ''  # 修改为原始txt文件所在的目录
mat_dir_eta = 'pixel_data/eta'       # 修改为eta.npy文件所在的目录
mat_dir_u = 'pixel_data/u'           # 修改为u.npy文件所在的目录
mat_dir_v = 'pixel_data/v'           # 修改为v.npy文件所在的目录

# 定义文件数量和数据尺寸
num_files = 1000
rows = 4096
cols = 2048

# 采样次数
num_samples = 100

# 初始化随机数种子以确保可重复性
np.random.seed(42)

# 随机生成采样索引
sample_files = np.random.randint(0, num_files, num_samples)
sample_rows = np.random.randint(1, rows + 1, num_samples)
sample_cols = np.random.randint(1, cols + 1, num_samples)

def test_variable(args):
    var_name, file_index, row_index, col_index = args
    print(f"Testing: var_name={var_name}, file_index={file_index}, row_index={row_index}, col_index={col_index}")
    error_count = 0
    error_count = 0
    original_file = os.path.join(original_dir, f"{var_name}_{file_index:05d}")
    mat_file = os.path.join(f'pixel_data/{var_name}', f"{var_name}_{file_index:05d}.npy")

    # 从txt文件中读取相应的值
    data_txt = np.loadtxt(original_file).reshape(rows, cols)
    value_txt = data_txt[row_index - 1, col_index - 1]  # 转换成0-index

    # 从.npy文件中读取相应的值
    data_mat = np.load(mat_file)
    value_mat = data_mat[row_index - 1, col_index - 1]  # 转换成0-index

    # 检查是否一致
    if value_txt != value_mat:
        print(f'Mismatch in file {var_name}_{file_index:05d}.npy at row {row_index}, col {col_index}: txt value = {value_txt}, mat value = {value_mat}')
        error_count = 1
    
    return error_count

def main():
    with ProcessPoolExecutor() as executor:
        tasks = []
        for var_name in ['eta', 'u', 'v']:
            for sample_file, sample_row, sample_col in zip(sample_files, sample_rows, sample_cols):
                tasks.append((var_name, sample_file, sample_row, sample_col))
        
        results = list(executor.map(test_variable, tasks))
        error_count = sum(results)
        print(f'Total mismatches: {error_count}')

if __name__ == "__main__":
    main()


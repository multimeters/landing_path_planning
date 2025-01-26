import os
import numpy as np
from scipy.io import loadmat, savemat
import math

def calculate_vector_angle(u, v):
    # 计算每个像素点的向量角度
    angles = np.arctan2(v, u)  # 计算角度
    return angles

def process_files(u_folder, v_folder, a_folder, start_idx, end_idx):
    for i in range(start_idx, end_idx + 1):
        u_filename = f"u_{i * 256}-{(i + 1) * 256 - 1}.mat"
        v_filename = f"v_{i * 256}-{(i + 1) * 256 - 1}.mat"
        a_filename = f"a_{i * 256}-{(i + 1) * 256 - 1}.mat"
        
        u_path = os.path.join(u_folder, u_filename)
        v_path = os.path.join(v_folder, v_filename)
        a_path = os.path.join(a_folder, a_filename)

        # 加载u和v数据
        u_data = loadmat(u_path)['data']  # 假设mat文件中的变量名为'u'
        v_data = loadmat(v_path)['data']  # 假设mat文件中的变量名为'v'
        
        # 检查数据形状是否匹配
        if u_data.shape != v_data.shape:
            raise ValueError(f"Data shape mismatch in files {u_filename} and {v_filename}")
        
        # 计算向量角度
        angles = calculate_vector_angle(u_data, v_data)
        
        # 保存角度数据
        savemat(a_path, {'data': angles})

        print(f"Processed and saved: {a_filename}")

if __name__ == "__main__":
    u_folder = '../data/processing/u_output'
    v_folder = '../data/processing/v_output'
    a_folder = '../data/processing/a_output'
    start_idx = 0
    end_idx = 0  # 根据文件范围设置

    if not os.path.exists(a_folder):
        os.makedirs(a_folder)

    process_files(u_folder, v_folder, a_folder, start_idx, end_idx)


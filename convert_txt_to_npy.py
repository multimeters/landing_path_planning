import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor

def convert_txt_to_mat(prefix, index):
    filename = f"{prefix}_{index:05d}"
    
    # 读取数据，假设文件为以空格分隔的文本文件
    data = np.loadtxt(filename)
    
    # 创建保存路径
    save_path = os.path.join('pixel_data', prefix)
    
    # 如果文件夹不存在，则创建
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    print(f"{prefix}_{index:05d}")
    # 保存为.npy文件
    np.save(os.path.join(save_path, f"{prefix}_{index:05d}.npy"), data)

def main():
    # 删除现有的并行池并启动新的并行池
    executor = ProcessPoolExecutor()
    
    # 并行执行函数convert_txt_to_mat
    for prefix in ['eta', 'u', 'v']:
        for index in range(1000):       
            executor.submit(convert_txt_to_mat, prefix, index)
    
    # 关闭并行池
    executor.shutdown()

if __name__ == "__main__":
    main()


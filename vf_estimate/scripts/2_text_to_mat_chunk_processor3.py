import numpy as np
import pandas as pd
import os
from multiprocessing import Pool, cpu_count
import gc
import scipy.io
from functools import partial

def read_chunk_from_file(input_folder, prefix, file_idx, chunk_idx, chunk_size=256, 
                         num_rows=4096, num_cols=2048):
    """
    从单个文件中读取指定 chunk（行范围为 [chunk_idx*chunk_size, (chunk_idx+1)*chunk_size)）。
    返回 shape = (chunk_size, num_cols) 的 numpy 数组。
    """
    filename = f"{prefix}_{file_idx:05d}"  # 假设文件名类似 eta_00001 这样的格式
    filepath = os.path.join(input_folder, filename)
    
    # 计算需要跳过的行数
    skip_rows = chunk_idx * chunk_size
    
    # 如果最后一个 chunk 可能不满 chunk_size，则需要额外判断
    # 但如果 4096 / chunk_size 是整数，就不需要。
    # 这里假设 4096 能被 chunk_size 整除。
    
    try:
        # 只读指定行
        # 注意：pandas 的 skiprows 可以是 int 或者 callable/list
        #       这里传 int 就表示跳过 skip_rows 行
        data_df = pd.read_csv(
            filepath,
            sep='\s+',
            header=None,
            dtype=np.float32,
            skiprows=skip_rows,
            nrows=chunk_size,
            engine='c',
            memory_map=True
        )
        
        # 验证数据形状
        if data_df.shape != (chunk_size, num_cols):
            print(f"文件 {filename} 的第 {chunk_idx} 个块形状不正确: {data_df.shape}")
            return None
        
        return data_df.values
    
    except Exception as e:
        print(f"读取文件 {filename} 的第 {chunk_idx} 个块时出错: {e}")
        return None

def combine_all_files_by_chunk(input_folder, output_folder, prefix, 
                               total_files=1024, chunk_size=256, 
                               num_rows=4096, num_cols=2048):
    """
    通过「块级并行处理」的方式，逐个块地读取所有文件并保存。
    """
    # 计算总块数
    num_chunks = num_rows // chunk_size
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 准备一个进程池（可根据需要在外层统一创建，或在函数内创建）
    pool = Pool(processes=cpu_count())
    
    for chunk_idx in range(num_chunks):
        print(f"==> 开始处理第 {chunk_idx+1}/{num_chunks} 个块...")
        
        # 使用偏函数固定住公共参数
        partial_func = partial(
            read_chunk_from_file, 
            input_folder, prefix, 
            chunk_idx=chunk_idx, 
            chunk_size=chunk_size, 
            num_rows=num_rows, 
            num_cols=num_cols
        )
        
        # 并行读取该 chunk 的所有文件
        file_indices = range(total_files)
        results = pool.map(partial_func, file_indices)
        
        # 先计算成功读取文件的数量
        valid_results = [res for res in results if res is not None]
        if len(valid_results) < total_files:
            print(f"警告：有 {total_files - len(valid_results)} 个文件读取失败，跳过这些文件。")
        
        # 如果全部都读取成功，合并为 (chunk_size, num_cols, total_files)
        # 如果有读取失败的，可以考虑是否要跳过或补0等
        combined_chunk = np.zeros((chunk_size, num_cols, total_files), dtype=np.float32)
        
        for file_idx, chunk_data in enumerate(results):
            if chunk_data is not None:
                combined_chunk[:, :, file_idx] = chunk_data
        
        # 保存该 chunk 到 .mat 文件
        start_row = chunk_idx * chunk_size
        end_row = start_row + chunk_size - 1
        output_filename = os.path.join(output_folder, f"{prefix}_{start_row}-{end_row}.mat")
        scipy.io.savemat(output_filename, {'data': combined_chunk})
        print(f"已保存 {output_filename}")
        
        # 释放内存
        del results
        del valid_results
        del combined_chunk
        gc.collect()
    
    # 关闭进程池
    pool.close()
    pool.join()
    print("全部块处理完成")

def main():
    # 文件夹路径和前缀
    folders = [
        '../data/source/output/',
        '../data/source/output/',
        '../data/source/output/'
    ]
    output_folders = ['../data/processing/e_output', '../data/processing/u_output', '../data/processing/v_output']
    prefixes = ['eta', 'u', 'v']
    
    total_files = 1024  # 如果不同，请更新
    chunk_size = 256
    num_rows = 256
    num_cols = 512

    for input_folder, output_folder, prefix in zip(folders, output_folders, prefixes):
        print(f"开始处理文件夹: {input_folder} 对应前缀: {prefix}")
        combine_all_files_by_chunk(
            input_folder=input_folder,
            output_folder=output_folder,
            prefix=prefix,
            total_files=total_files,
            chunk_size=chunk_size,
            num_rows=num_rows,
            num_cols=num_cols
        )
        print(f"完成处理文件夹: {input_folder}\n")

if __name__ == "__main__":
    main()

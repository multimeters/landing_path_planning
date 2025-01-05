import numpy as np
import scipy.io
import sys

def main():
    # 定义输入MAT文件和输出TXT文件的路径
    cvar_file = 'data/cvar_results_combined.mat'            # 向量长度文件
    directions_file = 'data/main_directions.mat'            # 向量方向文件
    
    u_output = 'data/combined_u_vector_field_k0.03.txt'      # 输出的U分量文件
    v_output = 'data/combined_v_vector_field_k0.03.txt'      # 输出的V分量文件
    
    # 加载cvar_results_combined.mat文件
    try:
        print(f"正在加载文件: {cvar_file} ...")
        cvar_data = scipy.io.loadmat(cvar_file)
        if 'cvar_image' not in cvar_data:
            print(f"错误: 文件 {cvar_file} 中未找到变量 'cvar_image'")
            sys.exit(1)
        cvar_image = cvar_data['cvar_image']           # 获取变量名为'cvar_image'的数据
        print(f"{cvar_file} 加载成功，数据维度: {cvar_image.shape}")
    except Exception as e:
        print(f"加载文件 {cvar_file} 时出错: {e}")
        sys.exit(1)
    
    # 加载main_directions.mat文件
    try:
        print(f"正在加载文件: {directions_file} ...")
        directions_data = scipy.io.loadmat(directions_file)
        if 'main_directions' not in directions_data:
            print(f"错误: 文件 {directions_file} 中未找到变量 'main_directions'")
            sys.exit(1)
        main_directions = directions_data['main_directions']   # 获取变量名为'main_directions'的数据
        print(f"{directions_file} 加载成功，数据维度: {main_directions.shape}")
    except Exception as e:
        print(f"加载文件 {directions_file} 时出错: {e}")
        sys.exit(1)
    
    # 检查矩阵维度是否为4096x2048
    expected_shape = (4096, 2048)
    if cvar_image.shape != expected_shape:
        print(f"错误: {cvar_file} 的维度为 {cvar_image.shape}，预期为 {expected_shape}")
        sys.exit(1)
    
    if main_directions.shape != expected_shape:
        print(f"错误: {directions_file} 的维度为 {main_directions.shape}，预期为 {expected_shape}")
        sys.exit(1)
    
    # 处理NaN值：如果cvar_image或main_directions中的值为NaN，则对应的u和v设为向右的单位向量 (1, 0)
    print("正在处理NaN值...")
    # 创建一个布尔掩码，标记cvar_image或main_directions中为NaN的位置
    nan_mask = np.isnan(cvar_image) | np.isnan(main_directions)
    
    # 将NaN位置的cvar_image和main_directions设为1和0，以便后续计算u和v
    cvar_image_clean = np.copy(cvar_image)
    main_directions_clean = np.copy(main_directions)
    cvar_image_clean[nan_mask] = 1.0
    main_directions_clean[nan_mask] = 0.0
    print("NaN值处理完成。")
    
    # 计算U和V分量
    print("正在计算U和V分量...")
    u = cvar_image_clean * np.cos(main_directions_clean)          # U分量 = 长度 * cos(方向)
    v = cvar_image_clean * np.sin(main_directions_clean)          # V分量 = 长度 * sin(方向)
    
    # 修改：将NaN位置的U和V分量设为向右的单位向量 (1, 0)
    u[nan_mask] = 1.0
    v[nan_mask] = 0.0
    print("U和V分量计算完成。")
    
    # 保存U分量到文本文件
    try:
        print(f"正在保存U分量到文件: {u_output} ...")
        np.savetxt(u_output, u, fmt='%.6f')            # 保存为浮点数，保留6位小数
        print(f"U分量已成功保存到 {u_output}")
    except Exception as e:
        print(f"保存U分量到文件 {u_output} 时出错: {e}")
        sys.exit(1)
    
    # 保存V分量到文本文件
    try:
        print(f"正在保存V分量到文件: {v_output} ...")
        np.savetxt(v_output, v, fmt='%.6f')            # 保存为浮点数，保留6位小数
        print(f"V分量已成功保存到 {v_output}")
    except Exception as e:
        print(f"保存V分量到文件 {v_output} 时出错: {e}")
        sys.exit(1)
    
    print("所有文件已成功转换和保存。")

if __name__ == '__main__':
    main()

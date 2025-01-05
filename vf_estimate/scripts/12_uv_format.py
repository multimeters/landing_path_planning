import scipy.io
import numpy as np
import os

def load_mat_file(filename, variable_name):
    """
    加载指定的 .mat 文件并提取指定的变量。
    
    :param filename: .mat 文件的路径
    :param variable_name: 要提取的变量名
    :return: 变量的 numpy 数组
    """
    try:
        mat = scipy.io.loadmat(filename)
        if variable_name in mat:
            return mat[variable_name]
        else:
            raise KeyError(f"变量 '{variable_name}' 不存在于文件 '{filename}' 中。")
    except FileNotFoundError:
        raise FileNotFoundError(f"文件 '{filename}' 未找到。")
    except Exception as e:
        raise e

def compute_uv(cvar_image, main_directions, direction_in_degrees=False):
    """
    根据长度和方向计算向量的 u 和 v 分量。
    如果 cvar_image 或 main_directions 中存在 NaN，则对应的 u 和 v 分量设为 0。
    
    :param cvar_image: 向量长度的二维数组
    :param main_directions: 向量方向的二维数组
    :param direction_in_degrees: 如果方向是以度为单位，则设置为 True
    :return: u 和 v 分量的二维数组
    """
    # 创建掩码，标记 cvar_image 或 main_directions 中的 NaN 值
    nan_mask = np.isnan(cvar_image) | np.isnan(main_directions)
    
    if direction_in_degrees:
        # 将角度转换为弧度
        main_directions = np.deg2rad(main_directions)
    
    # 计算 u 和 v 分量
    u = cvar_image * np.cos(main_directions)
    v = cvar_image * np.sin(main_directions)
    
    # 将 NaN 位置的 u 和 v 设为 0
    u[nan_mask] = 0
    v[nan_mask] = 0
    
    return u, v

def save_to_txt(data, filename):
    """
    将二维数组保存为文本文件，每行对应数组的一行，元素之间用空格分隔。
    
    :param data: 要保存的二维数组
    :param filename: 输出文件名（包含路径）
    """
    # 确保保存目录存在
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    np.savetxt(filename, data, fmt='%.6f')  # 保留6位小数，根据需要调整

def main():
    # 文件名和变量名
    cvar_mat_file = '../data/output/cvar_results_combined.mat'
    cvar_variable = 'cvar_image'
    
    directions_mat_file = '../data/output/main_directions.mat'
    directions_variable = 'main_directions'
    
    # 加载数据
    print("加载 cvar_image 数据...")
    cvar_image = load_mat_file(cvar_mat_file, cvar_variable)
    
    print("加载 main_directions 数据...")
    main_directions = load_mat_file(directions_mat_file, directions_variable)
    
    # 验证维度
    expected_shape = (256, 512)
    if cvar_image.shape != expected_shape:
        raise ValueError(f"cvar_image 的维度为 {cvar_image.shape}，预期为 {expected_shape}。")
    if main_directions.shape != expected_shape:
        raise ValueError(f"main_directions 的维度为 {main_directions.shape}，预期为 {expected_shape}。")
    
    print("计算 u 和 v 分量（处理 NaN 值）...")
    # 如果方向是以度为单位，请将 direction_in_degrees 设置为 True
    u, v = compute_uv(cvar_image, main_directions, direction_in_degrees=False)
    
    # 定义输出目录和文件路径
    output_dir = os.path.join('../data/output/', 'vf')
    u_output_file = os.path.join(output_dir, 'u_vf.txt')
    v_output_file = os.path.join(output_dir, 'v_vf.txt')
    
    # 保存结果
    print(f"保存 u 分量到 {u_output_file} ...")
    save_to_txt(u, u_output_file)
    
    print(f"保存 v 分量到 {v_output_file} ...")
    save_to_txt(v, v_output_file)
    
    print("完成！")

if __name__ == "__main__":
    main()

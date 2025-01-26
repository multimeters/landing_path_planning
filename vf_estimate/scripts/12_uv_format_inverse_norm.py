import scipy.io
import numpy as np
import os
import matplotlib

# 设置 matplotlib 使用 'Agg' 后端，以避免 Qt 平台插件错误
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

def compute_uv(cvar_image, main_directions, invert=False, direction_in_degrees=False):
    """
    根据长度和方向计算向量的 u 和 v 分量。
    如果 invert 为 True，则使用 cvar_image 的倒数作为向量长度。
    如果 cvar_image 或 main_directions 中存在 NaN，则对应的 u 和 v 分量设为 0。
    
    :param cvar_image: 向量长度的二维数组
    :param main_directions: 向量方向的二维数组
    :param invert: 是否使用 cvar_image 的倒数作为向量长度
    :param direction_in_degrees: 如果方向是以度为单位，则设置为 True
    :return: u 和 v 分量的二维数组
    """
    # 创建掩码，标记 cvar_image 或 main_directions 中的 NaN 值
    nan_mask = np.isnan(cvar_image) | np.isnan(main_directions)
    
    if invert:
        epsilon = 1e-6  # 防止除零
        cvar_image = 1.0 / (cvar_image + epsilon)
    
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

def plot_magnitude_distribution(u, v, output_dir, magnitude_threshold=0.1):
    """
    计算向量的幅值，并绘制其分布图（直方图），排除幅值小于指定阈值的向量。
    
    :param u: 向量的 u 分量
    :param v: 向量的 v 分量
    :param output_dir: 图像的保存目录
    :param magnitude_threshold: 幅值阈值，小于该值的向量将被排除
    :return: 幅值数组（过滤后的）
    """
    # 计算幅值
    magnitude = np.sqrt(u**2 + v**2)
    
    # 过滤掉幅值小于阈值和 NaN 值
    filtered_magnitude = magnitude[~np.isnan(magnitude) & (magnitude >= magnitude_threshold)]
    
    # 绘制直方图
    plt.figure(figsize=(10, 6))
    plt.hist(filtered_magnitude.flatten(), bins=100, color='blue', edgecolor='black', alpha=0.7)
    plt.title(f'向量幅值分布 (幅值 ≥ {magnitude_threshold})')
    plt.xlabel('幅值')
    plt.ylabel('频数')
    plt.grid(True)
    
    # 确保保存目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存图像
    plot_filename = os.path.join(output_dir, 'magnitude_distribution.png')
    plt.savefig(plot_filename)
    plt.close()
    
    print(f"向量幅值分布图已保存到 {plot_filename}.")
    
    return filtered_magnitude

def normalize_magnitude(filtered_magnitude, method='log', log_epsilon=1e-6, quantile=0.95):
    """
    根据指定的方法对幅值进行归一化。
    
    :param filtered_magnitude: 过滤后的幅值数组
    :param method: 归一化方法，支持 'log' 和 'quantile'
    :param log_epsilon: 对数归一化时的平滑因子，防止取对数时出现负无穷
    :param quantile: 分位数归一化时使用的分位数（0 < quantile < 1）
    :return: 归一化因子（浮点数）
    """
    if method == 'log':
        # 对数归一化
        normalized_magnitude = np.log(filtered_magnitude + log_epsilon)
        # 选择一个合理的归一化因子，例如对数均值
        normalization_factor = np.mean(normalized_magnitude)
        print(f"使用对数归一化，归一化因子为对数均值: {normalization_factor:.4f}")
    elif method == 'quantile':
        # 分位数归一化
        normalization_factor = np.quantile(filtered_magnitude, quantile)
        normalized_magnitude = filtered_magnitude / normalization_factor
        print(f"使用分位数归一化，归一化因子（{quantile*100}分位数）: {normalization_factor:.4f}")
    else:
        raise ValueError("Unsupported normalization method. Choose 'log' or 'quantile'.")
    
    return normalization_factor

def apply_normalization(u, v, normalization_factor, method='log', log_epsilon=1e-6):
    """
    根据归一化因子对 u 和 v 分量进行归一化。
    
    :param u: 向量的 u 分量
    :param v: 向量的 v 分量
    :param normalization_factor: 归一化因子
    :param method: 归一化方法，支持 'log' 和 'quantile'
    :param log_epsilon: 对数归一化时的平滑因子，防止取对数时出现负无穷
    :return: 归一化后的 u 和 v 分量
    """
    magnitude = np.sqrt(u**2 + v**2)
    mask = (~np.isnan(magnitude)) & (magnitude >= 0.1)
    
    if method == 'log':
        # 对数归一化
        normalized_magnitude = np.log(magnitude[mask] + log_epsilon) / normalization_factor
    elif method == 'quantile':
        # 分位数归一化
        normalized_magnitude = magnitude[mask] / normalization_factor
    else:
        raise ValueError("Unsupported normalization method. Choose 'log' or 'quantile'.")
    
    # 计算归一化后的 u 和 v
    u_normalized = np.zeros_like(u)
    v_normalized = np.zeros_like(v)
    u_normalized[mask] = (u[mask] / magnitude[mask]) * normalized_magnitude
    v_normalized[mask] = (v[mask] / magnitude[mask]) * normalized_magnitude
    
    return u_normalized, v_normalized

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
    
    print("计算 u 和 v 分量前，统计向量幅值分布...")
    # 如果方向是以度为单位，请将 invert 设置为 True
    # 根据您的需求，决定是否取倒数
    invert = False  # 设置为 True 以使用 cvar_image 的倒数
    u_temp, v_temp = compute_uv(cvar_image, main_directions, invert=invert, direction_in_degrees=False)
    
    # 定义输出目录
    output_dir = os.path.join('../data/output/', 'vf')
    
    # 设置幅值阈值
    magnitude_threshold = 0.1  # 幅值小于0.1的向量将被排除
    
    # 绘制幅值分布图，排除幅值小于阈值的向量
    filtered_magnitude = plot_magnitude_distribution(u_temp, v_temp, output_dir, magnitude_threshold=magnitude_threshold)
    
    # 选择归一化方法
    normalization_method = 'quantile'  # 选择 'log' 或 'quantile'
    
    # 计算归一化因子
    normalization_factor = normalize_magnitude(filtered_magnitude, method=normalization_method)
    
    print("应用归一化到 u 和 v 分量...")
    # 归一化 u 和 v
    u_normalized, v_normalized = apply_normalization(u_temp, v_temp, normalization_factor, method=normalization_method)
    
    # 定义输出文件路径
    u_output_file = os.path.join(output_dir, 'u_vf_normalized.txt')
    v_output_file = os.path.join(output_dir, 'v_vf_normalized.txt')
    
    # 保存归一化后的结果
    print(f"保存归一化后的 u 分量到 {u_output_file} ...")
    save_to_txt(u_normalized, u_output_file)
    
    print(f"保存归一化后的 v 分量到 {v_output_file} ...")
    save_to_txt(v_normalized, v_output_file)
    
    # 可选：再次绘制归一化后的幅值分布
    plot_magnitude_distribution(u_normalized, v_normalized, output_dir, magnitude_threshold=magnitude_threshold)
    
    print("完成！")

if __name__ == "__main__":
    main()

import numpy as np
import scipy.io as sio
from scipy.stats import norm
import os
import glob
import re

import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt

def load_mat_file(file_path):
    """
    加载.mat文件并提取相关变量。
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件 {file_path} 不存在。")
    
    try:
        mat = sio.loadmat(file_path)
    except NotImplementedError:
        # 如果是基于HDF5的.mat文件，使用h5py加载
        import h5py
        with h5py.File(file_path, 'r') as f:
            mat = {k: np.array(v) for k, v in f.items()}
    
    required_vars = ['var_x3', 'var_x4', 'var_x5', 'cov_x3_x4', 'cov_x3_x5', 'cov_x4_x5']
    for var in required_vars:
        if var not in mat:
            raise KeyError(f"变量 '{var}' 不存在于 {file_path} 中。")
    
    var_x3 = mat['var_x3']
    var_x4 = mat['var_x4']
    var_x5 = mat['var_x5']
    cov_x3_x4 = mat['cov_x3_x4']
    cov_x3_x5 = mat['cov_x3_x5']
    cov_x4_x5 = mat['cov_x4_x5']
    
    return var_x3, var_x4, var_x5, cov_x3_x4, cov_x3_x5, cov_x4_x5

def construct_covariance_matrices(var_x3, var_x4, var_x5, cov_x3_x4, cov_x3_x5, cov_x4_x5):
    """
    构建每个像素的3x3协方差矩阵。
    返回形状为 (num_pixels, 3, 3) 的numpy数组。
    """
    num_rows, num_cols = var_x3.shape
    num_pixels = num_rows * num_cols
    
    # 初始化协方差矩阵数组
    cov_matrices = np.zeros((num_pixels, 3, 3))
    
    # 填充协方差矩阵
    cov_matrices[:, 0, 0] = var_x3.flatten()
    cov_matrices[:, 1, 1] = var_x4.flatten()
    cov_matrices[:, 2, 2] = var_x5.flatten()
    cov_matrices[:, 0, 1] = cov_x3_x4.flatten()
    cov_matrices[:, 1, 0] = cov_x3_x4.flatten()
    cov_matrices[:, 0, 2] = cov_x3_x5.flatten()
    cov_matrices[:, 2, 0] = cov_x3_x5.flatten()
    cov_matrices[:, 1, 2] = cov_x4_x5.flatten()
    cov_matrices[:, 2, 1] = cov_x4_x5.flatten()
    
    return cov_matrices, num_rows, num_cols

def compute_cvar(cov_matrices, w, alpha=0.95):
    """
    计算每个像素的CVaR。
    
    参数：
    - cov_matrices: 形状为 (num_pixels, 3, 3) 的协方差矩阵数组
    - w: 权重向量，形状为 (3,)
    - alpha: 置信水平，默认0.95
    
    返回：
    - cvar_values: 形状为 (num_pixels,) 的CVaR值
    """
    # 计算 sigma_L = w^T * Σ * w
    sigma_L_all = np.sum(cov_matrices * np.outer(w, w), axis=(1,2))
    
    # 计算CVaR
    z_alpha = norm.ppf(alpha)
    cvar_values = (norm.pdf(z_alpha) / (1 - alpha)) * np.sqrt(sigma_L_all)
    
    return cvar_values

def visualize_cvar(cvar_image, title='CVaR per Pixel', cmap='hot'):
    """
    可视化CVaR图像。
    """
    plt.figure(figsize=(20, 15))  # 增大图像尺寸以适应拼接后的大图
    plt.imshow(cvar_image, cmap=cmap)
    plt.colorbar(label='CVaR')
    plt.title(title)
    plt.xlabel('Columns')
    plt.ylabel('Rows')
    plt.show()

def save_cvar_image(cvar_image, output_path, cmap='hot'):
    """
    保存CVaR图像为PNG文件。
    """
    plt.figure(figsize=(20, 15))  # 增大图像尺寸以适应拼接后的大图
    plt.imshow(cvar_image, cmap=cmap)
    plt.colorbar(label='CVaR')
    plt.title('Combined CVaR per Pixel')
    plt.xlabel('Columns')
    plt.ylabel('Rows')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"CVaR图像已保存到 {output_path}")

def save_cvar_mat(cvar_image, output_mat_path):
    """
    保存CVaR图像到.mat文件。
    """
    sio.savemat(output_mat_path, {'cvar_image': cvar_image})
    print(f"CVaR结果已保存到 {output_mat_path}")

def extract_row_range(file_name):
    """
    从文件名中提取行范围。
    假设文件名格式为 'variance_covariance_results_per_pixel-<start>-<end>.mat'
    """
    pattern = r'variance_covariance_results_per_pixel-(\d+)-(\d+)\.mat'
    match = re.search(pattern, file_name)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return start, end
    else:
        raise ValueError(f"文件名 {file_name} 不符合预期格式。")

def arrange_images_grid(images, grid_size=(4,4)):
    """
    将图像列表按指定网格大小拼接。
    
    参数：
    - images: list of 2D numpy arrays
    - grid_size: tuple, e.g., (4,4) for 4 rows and 4 columns
    
    返回：
    - concatenated_image: 2D numpy array
    """
    if len(images) != grid_size[0] * grid_size[1]:
        raise ValueError(f"需要 {grid_size[0] * grid_size[1]} 个图像，但提供了 {len(images)} 个。")
    
    # 按行分组
    rows = []
    for i in range(grid_size[0]):
        row_images = images[i * grid_size[1] : (i + 1) * grid_size[1]]
        row = np.hstack(row_images)
        rows.append(row)
    
    # 垂直拼接所有行
    concatenated_image = np.vstack(rows)
    
    return concatenated_image

def process_files(input_pattern, w, alpha=0.95, visualize=False, grid_size=(4,4)):
    """
    处理多个.mat文件，计算CVaR，拼接结果，并保存最终的拼接图像和.mat文件。
    
    参数：
    - input_pattern: 字符串，匹配输入.mat文件的模式
    - w: numpy数组，权重向量
    - alpha: float，置信水平
    - visualize: bool，是否可视化每个CVaR图像
    - grid_size: tuple, 拼接网格大小 (rows, cols)
    """
    # 找到所有匹配的.mat文件，并按行范围排序
    input_files = sorted(glob.glob(input_pattern), key=lambda x: extract_row_range(os.path.basename(x))[0])
    
    expected_num_files = grid_size[0] * grid_size[1]
    if len(input_files) != expected_num_files:
        print(f"警告: 预期处理 {expected_num_files} 个文件，但找到 {len(input_files)} 个文件。")
    
    print(f"找到 {len(input_files)} 个文件需要处理。")
    
    cvar_images = []  # 存储每个文件的CVaR图像
    
    for input_mat_file in input_files:
        try:
            print(f"\n正在处理文件: {input_mat_file}")
            
            # 加载数据
            print("加载.mat文件...")
            var_x3, var_x4, var_x5, cov_x3_x4, cov_x3_x5, cov_x4_x5 = load_mat_file(input_mat_file)
            
            # 构建协方差矩阵
            print("构建协方差矩阵...")
            cov_matrices, num_rows, num_cols = construct_covariance_matrices(var_x3, var_x4, var_x5, cov_x3_x4, cov_x3_x5, cov_x4_x5)
            
            # 计算CVaR
            print("计算CVaR...")
            cvar_values = compute_cvar(cov_matrices, w, alpha)
            
            # 将CVaR值重新调整为图像形状
            cvar_image = cvar_values.reshape(num_rows, num_cols)
            cvar_images.append(cvar_image)
            
            # 可视化
            if visualize:
                print("可视化CVaR图像...")
                visualize_cvar(cvar_image, title=f'CVaR per Pixel - {os.path.basename(input_mat_file)}')
            
            print("文件处理完成。")
            
        except Exception as e:
            print(f"处理文件 {input_mat_file} 时出错: {e}")
    
    # 检查是否有足够的图像进行拼接
    if len(cvar_images) < expected_num_files:
        print(f"注意: 只有 {len(cvar_images)} 个文件成功处理，无法完全拼接 {expected_num_files} 个文件的图像。")
    
    if not cvar_images:
        print("没有成功处理任何文件，无法生成拼接图像。")
        return
    
    # 按行范围顺序纵向拼接所有CVaR图像
    print("\n拼接所有CVaR图像（纵向）...")
    try:
        concatenated_cvar_image = np.vstack(cvar_images)
        print(f"拼接后图像形状: {concatenated_cvar_image.shape}")
    except Exception as e:
        print(f"拼接图像时出错: {e}")
        return
    
    # 保存拼接后的图像
    output_image_path = '../data/output/cvar_image_combined.png'
    print(f"保存拼接后的CVaR图像到 {output_image_path}...")
    save_cvar_image(concatenated_cvar_image, output_image_path)
    
    # 保存拼接后的.mat文件
    output_mat_file = '../data/output/cvar_results_combined.mat'
    print(f"保存拼接后的CVaR结果到 {output_mat_file}...")
    save_cvar_mat(cvar_image, output_mat_file)
    
    print("\n所有文件已成功处理和保存。")

def main():
    # 配置参数
    input_folder = '../data/processing/conv_cal'  # 输入文件夹
    input_pattern = os.path.join(input_folder, 'variance_covariance_results_per_pixel-*.mat')  # 输入文件模式
    alpha = 0.95  # 置信水平
    visualize = False  # 是否可视化每个CVaR图像
    grid_size = (4,4)  # 4行4列的网格
    
    # 定义权重向量
    # 请根据实际情况设置w1, w2, w3
    # 例如，均匀权重:
    w = np.array([1.0, 1.0, 1.0])
    
    # 或者根据具体需求设置不同权重:
    # w = np.array([w1, w2, w3])
    
    # 检查输入文件夹是否存在
    if not os.path.isdir(input_folder):
        print(f"输入文件夹 '{input_folder}' 不存在。请确保文件夹存在并包含所需的.mat文件。")
        return
    
    # 处理所有文件并拼接结果
    process_files(input_pattern, w, alpha, visualize, grid_size)

if __name__ == "__main__":
    main()

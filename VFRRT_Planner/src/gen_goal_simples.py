import numpy as np
from skimage import measure
import os

def read_depth_data(file_path, rows=256, cols=512):
    """
    读取深度数据文件，并返回一个二维numpy数组。
    假设数据以空格或逗号分隔。
    """
    try:
        depth = np.loadtxt(file_path, delimiter=None)  # 自动检测分隔符
        if depth.shape != (rows, cols):
            raise ValueError(f"数据维度与预期不符: {depth.shape} vs {(rows, cols)}")
        return depth
    except Exception as e:
        print(f"读取深度数据失败: {e}")
        raise

def extract_contour(depth, level=0.0):
    """
    提取指定深度的等值线。
    返回轮廓的x和y坐标数组。
    """
    contours = measure.find_contours(depth, level=level)
    if not contours:
        raise ValueError(f"未找到深度为{level}的等值线。")
    # 如果有多个轮廓，可以选择最长的一个作为主要海岸线
    contour = max(contours, key=lambda x: len(x))
    return contour

def convert_to_real_coordinates(contour, x_resolution=1.0, y_resolution=2.0):
    """
    将像素坐标转换为实际坐标。
    x: 列索引 * x_resolution
    y: 行索引 * y_resolution
    """
    x = contour[:, 1] * x_resolution
    y = contour[:, 0] * y_resolution
    return x, y

def sample_along_contour(x, y, sample_distance=5.0):
    """
    沿着轮廓路径每隔sample_distance米采样一个点。
    返回采样点的x和y坐标列表。
    """
    points = np.vstack((x, y)).T
    # 计算每两点之间的距离
    deltas = np.diff(points, axis=0)
    seg_lengths = np.hypot(deltas[:,0], deltas[:,1])
    cumulative_lengths = np.insert(np.cumsum(seg_lengths), 0, 0)

    total_length = cumulative_lengths[-1]
    num_samples = int(total_length // sample_distance) + 1
    sample_points = []

    for i in range(num_samples):
        target_length = i * sample_distance
        if target_length > total_length:
            break
        # 找到目标长度所在的段
        idx = np.searchsorted(cumulative_lengths, target_length) - 1
        if idx >= len(seg_lengths):
            idx = len(seg_lengths) - 1
        # 计算在该段的比例
        segment_start = cumulative_lengths[idx]
        segment_end = cumulative_lengths[idx + 1]
        segment_ratio = (target_length - segment_start) / seg_lengths[idx]
        # 线性插值计算采样点
        sample_x = points[idx,0] + deltas[idx,0] * segment_ratio
        sample_y = points[idx,1] + deltas[idx,1] * segment_ratio
        sample_points.append((sample_x, sample_y))

    return sample_points

def save_sampled_points(sample_points, output_file):
    """
    将采样点保存到文本文件中，每行包含x和y坐标，逗号分隔。
    """
    try:
        with open(output_file, 'w') as f:
            for point in sample_points:
                f.write(f"{point[0]:.2f} {point[1]:.2f}\n")
        print(f"采样点已保存到 {output_file}")
    except Exception as e:
        print(f"保存采样点失败: {e}")
        raise

def main():
    depth_file = "../data/dep_shoal_inlet_subset_swapped.txt"
    output_file = "../data/sampled_coastline.txt"

    if not os.path.exists(depth_file):
        print(f"深度数据文件 {depth_file} 不存在。")
        return

    # 步骤1：读取深度数据
    depth = read_depth_data(depth_file)

    # 步骤2：提取深度为0.75的等值线
    contour = extract_contour(depth, level=0.75)

    # 步骤3：转换为实际坐标
    x, y = convert_to_real_coordinates(contour)

    # 步骤4：沿海岸线每隔5米采样
    sample_points = sample_along_contour(x, y, sample_distance=5.0)

    # 步骤5：保存采样点到文本文件
    save_sampled_points(sample_points, output_file)

if __name__ == "__main__":
    main()

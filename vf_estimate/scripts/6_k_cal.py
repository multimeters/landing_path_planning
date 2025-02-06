import numpy as np
from scipy.io import savemat

def compute_expression(input_txt, output_mat):
    # 定义常数
    U_kmh = 10  # 速度，单位：公里/小时
    U = U_kmh * 1000 / 3600  # 转换为米/秒
    mu = 0  # 角度，单位：弧度
    g = 9.8  # 重力加速度，单位：米/秒²

    # 加载深度数据
    try:
        # 假设数据以空格、逗号或制表符分隔
        h = np.loadtxt(input_txt, delimiter=None)  # 自动检测分隔符
    except Exception as e:
        print(f"加载文件 {input_txt} 时出错: {e}")
        return

    # 确保h是一个二维数组
    if h.ndim != 2:
        print(f"输入数据不是二维数组，当前维度为 {h.ndim}")
        return

    # 计算cos(mu)，由于mu=0，cos(mu)=1
    cos_mu = np.cos(mu)

    # 计算sqrt(g * h)
    with np.errstate(invalid='ignore'):  # 忽略无效值的警告
        sqrt_gh = np.sqrt(g * h)

    # 计算表达式
    with np.errstate(divide='ignore', invalid='ignore'):
        result = 1 - (U * cos_mu) / sqrt_gh

    # 对于h < 0.75，设置结果为NaN
    result[h < 0.75] = np.nan

    # 保存结果到.mat文件
    try:
        savemat(output_mat, {'k': result})
        print(f"计算结果已保存到 {output_mat}")
    except Exception as e:
        print(f"保存文件 {output_mat} 时出错: {e}")

if __name__ == "__main__":
    # 输入和输出文件名
    input_txt = '../data/source/dep_shoal_inlet_subset_swapped.txt'       # 替换为您的输入文件名
    output_mat = '../data/processing/k.mat' # 替换为您希望的输出文件名

    compute_expression(input_txt, output_mat)

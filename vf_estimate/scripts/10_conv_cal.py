import numpy as np
import scipy.io as sio
import sys
import os

def main():
    try:
        print("=== 开始处理所有数据块 ===")

        # 定义阈值
        max3 = 3.0   # Heave 阈值，单位：米
        max4 = 0.52  # Roll 阈值（30度）
        max5 = 0.52  # Pitch 阈值（30度）
        print("阈值已定义。")

        # 定义每种文件类型的基础文件夹路径
        # 请根据实际情况修改以下路径
        heave_base_dir = '../data/processing/output_heave_rao'    # 替换为heave文件所在的文件夹路径
        pitch_base_dir = '../data/processing/output_pitch_rao'    # 替换为pitch文件所在的文件夹路径
        energy_base_dir = '../data/processing/fft_results'        # 替换为energy文件所在的文件夹路径
        mask_base_dir = '../data/processing/mask'                 # 替换为mask文件所在的文件夹路径
        output_base_dir = '../data/processing/conv_cal'           # 替换为结果文件保存的文件夹路径

        # 确保输出文件夹存在，如果不存在则创建
        if not os.path.isdir(output_base_dir):
            os.makedirs(output_base_dir)
            print(f"已创建输出文件夹: {output_base_dir}")

        # 定义数据块的数量和大小
        num_blocks = 1
        block_size = 256  # 每个块的大小为256（如0-255, 256-511, ...）

        for block in range(num_blocks):
            start = block * block_size
            end = (block + 1) * block_size - 1
            block_label = f"{start}-{end}"
            print(f"\n--- 处理数据块 {block + 1}/{num_blocks}: {block_label} ---")

            # 构建完整的文件路径
            heave_filename = os.path.join(heave_base_dir, f'heave_rao-{block_label}.mat')
            pitch_filename = os.path.join(pitch_base_dir, f'pitch_rao-{block_label}.mat')
            energy_filename = os.path.join(energy_base_dir, f'energy_block_{block_label}.mat')  # 修改文件扩展名为 .mat
            mask_filename = os.path.join(mask_base_dir, f'mask.mat')
            output_filename = os.path.join(output_base_dir, f'variance_covariance_results_per_pixel-{block_label}.mat')

            # 检查输入文件是否存在
            print("检查所需文件是否存在...")
            for fname in [heave_filename, pitch_filename, energy_filename, mask_filename]:
                if not os.path.isfile(fname):
                    raise FileNotFoundError(f"文件 '{fname}' 不存在。请检查路径和文件名。")
            print("所有所需文件均存在。")

            # 加载 RAO_3 从 heave_rao-<label>.mat
            print(f"从 '{heave_filename}' 加载 RAO_3...")
            heave_mat = sio.loadmat(heave_filename)
            RAO_3 = heave_mat.get('heave_rao')
            if RAO_3 is None:
                raise KeyError(f"在 '{heave_filename}' 中未找到变量 'heave_rao'。")
            print(f"RAO_3 加载成功，形状为 {RAO_3.shape}。")

            # 初始化 RAO_4 作为全零数组，形状与 RAO_3 相同
            print("初始化 RAO_4 为全零数组...")
            RAO_4 = np.zeros_like(RAO_3)
            print(f"RAO_4 初始化完成，形状为 {RAO_4.shape}。")

            # 加载 RAO_5 从 pitch_rao-<label>.mat
            print(f"从 '{pitch_filename}' 加载 RAO_5...")
            pitch_mat = sio.loadmat(pitch_filename)
            RAO_5 = pitch_mat.get('pitch_rao')
            if RAO_5 is None:
                raise KeyError(f"在 '{pitch_filename}' 中未找到变量 'pitch_rao'。")
            print(f"RAO_5 加载成功，形状为 {RAO_5.shape}。")

            # 加载 energy_density 从 energy_block-<label>.mat（修改为从 .mat 文件加载）
            print(f"从 '{energy_filename}' 加载 energy_density...")
            energy_mat = sio.loadmat(energy_filename)
            energy_density = energy_mat.get('energy_density_half')  # 根据实际变量名修改
            if energy_density is None:
                raise KeyError(f"在 '{energy_filename}' 中未找到变量 'energy_density_half'。请确认变量名。")
            print(f"energy_density 加载成功，形状为 {energy_density.shape}。")

            # 加载 mask 从 mask.mat
            print(f"从 '{mask_filename}' 加载掩膜...")
            mask_mat = sio.loadmat(mask_filename)
            mask = mask_mat.get('mask')
            if mask is None:
                raise KeyError(f"在 '{mask_filename}' 中未找到变量 'mask'。")
            print(f"掩膜加载成功，形状为 {mask.shape}。")

            # 确保掩膜为二值（0 或 1）
            print("将掩膜转换为布尔类型...")
            mask_bool = mask.astype(bool)  # 转换为布尔类型：1 为 True，0 为 False
            print("掩膜转换完成。")

            # 确保所有数组形状兼容
            print("验证所有数组的形状是否兼容...")
            if not (RAO_3.shape == RAO_4.shape == RAO_5.shape == energy_density.shape):
                raise ValueError("RAO 和 energy_density 数组的形状必须相同。")
            if mask_bool.shape != RAO_3.shape[:2]:
                raise ValueError("掩膜的形状必须与 RAO 和 energy_density 数组的前两维相同。")
            print("数组形状验证通过。")

            # 应用掩膜
            print("应用掩膜到数据...")

            # 使用掩膜将不需要处理的像素的数据设置为 0
            RAO_3_masked = np.where(mask_bool[:, :, np.newaxis], RAO_3, 0)
            RAO_4_masked = np.where(mask_bool[:, :, np.newaxis], RAO_4, 0)
            RAO_5_masked = np.where(mask_bool[:, :, np.newaxis], RAO_5, 0)
            energy_density_masked = np.where(mask_bool[:, :, np.newaxis], energy_density, 0)
            print("掩膜应用完成。")

            # 计算方差
            print("开始计算方差...")
            var_x3 = np.sum(energy_density_masked * (RAO_3_masked ** 2), axis=2) / (max3 ** 2)  # Heave 方差
            var_x4 = np.sum(energy_density_masked * (RAO_4_masked ** 2), axis=2) / (max4 ** 2)  # Roll 方差
            var_x5 = np.sum(energy_density_masked * (RAO_5_masked ** 2), axis=2) / (max5 ** 2)  # Pitch 方差
            print("方差计算完成。")

            # 计算协方差
            print("开始计算协方差...")
            cov_x3_x4 = np.sum(energy_density_masked * (RAO_3_masked * RAO_4_masked), axis=2) / (max3 * max4)  # Heave-Roll 协方差
            cov_x3_x5 = np.sum(energy_density_masked * (RAO_3_masked * RAO_5_masked), axis=2) / (max3 * max5)  # Heave-Pitch 协方差
            cov_x4_x5 = np.sum(energy_density_masked * (RAO_4_masked * RAO_5_masked), axis=2) / (max4 * max5)  # Roll-Pitch 协方差
            print("协方差计算完成。")

            # 将掩膜为 False 的像素的方差和协方差设置为 NaN
            print("将掩膜为 False 的像素的方差和协方差设置为 NaN...")
            var_x3 = var_x3.astype(np.float32)
            var_x4 = var_x4.astype(np.float32)
            var_x5 = var_x5.astype(np.float32)
            cov_x3_x4 = cov_x3_x4.astype(np.float32)
            cov_x3_x5 = cov_x3_x5.astype(np.float32)
            cov_x4_x5 = cov_x4_x5.astype(np.float32)

            # 使用 np.where 进行替换
            var_x3 = np.where(mask_bool, var_x3, np.nan)
            var_x4 = np.where(mask_bool, var_x4, np.nan)
            var_x5 = np.where(mask_bool, var_x5, np.nan)
            cov_x3_x4 = np.where(mask_bool, cov_x3_x4, np.nan)
            cov_x3_x5 = np.where(mask_bool, cov_x3_x5, np.nan)
            cov_x4_x5 = np.where(mask_bool, cov_x4_x5, np.nan)
            print("设置完成。")

            # 打印部分结果以供验证（可选）
            print("\n=== 部分计算结果 ===")
            print(f"var_x3（Heave 方差）示例值：{var_x3[0,0]}")
            print(f"var_x4（Roll 方差）示例值：{var_x4[0,0]}")
            print(f"var_x5（Pitch 方差）示例值：{var_x5[0,0]}")
            print(f"cov_x3_x4（Heave-Roll 协方差）示例值：{cov_x3_x4[0,0]}")
            print(f"cov_x3_x5（Heave-Pitch 协方差）示例值：{cov_x3_x5[0,0]}")
            print(f"cov_x4_x5（Roll-Pitch 协方差）示例值：{cov_x4_x5[0,0]}")

            # 准备结果字典
            print("\n准备将结果保存到 .mat 文件中...")
            results = {
                'var_x3': var_x3,
                'var_x4': var_x4,
                'var_x5': var_x5,
                'cov_x3_x4': cov_x3_x4,
                'cov_x3_x5': cov_x3_x5,
                'cov_x4_x5': cov_x4_x5
            }

            # 保存结果到 .mat 文件
            print(f"将结果保存到 '{output_filename}'...")
            sio.savemat(output_filename, results)
            print("结果保存成功。")

        print("\n=== 所有数据块处理完成 ===")

    except Exception as e:
        print(f"发生错误：{e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

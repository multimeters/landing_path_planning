import numpy as np
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool, cpu_count
from functools import partial

def generate_frame(i, nx, ny, skip, save_dir, dpi):
    """
    生成单个帧的函数。
    """
    try:
        # 拼接文件名
        eta_filename = f"/media/ganzhi/Crucial X6/output_waiyuan/output/eta_{i:05d}"
        u_filename   = f"/media/ganzhi/Crucial X6/output_waiyuan/output/u_{i:05d}"
        v_filename   = f"/media/ganzhi/Crucial X6/output_waiyuan/output/v_{i:05d}"
        
        # 读取数据
        eta = np.loadtxt(eta_filename).reshape(ny, nx)
        u   = np.loadtxt(u_filename).reshape(ny, nx)
        v   = np.loadtxt(v_filename).reshape(ny, nx)
        
        # 创建一个新的绘图对象，以确保线程安全
        fig, ax = plt.subplots(figsize=(nx / dpi, ny / dpi), dpi=dpi)
        
        # 绘制 eta 的灰度图
        im = ax.imshow(
            eta, 
            cmap='gray', 
            origin='lower', 
            extent=[0, nx, 0, ny]
        )
        
        # 绘制 (u, v) 的矢量场
        X = np.arange(0, nx, skip)
        Y = np.arange(0, ny, skip)
        xq, yq = np.meshgrid(X, Y)
        u_s = u[yq, xq]
        v_s = v[yq, xq]
        
        ax.quiver(
            xq, yq, 
            u_s, v_s, 
            color='red', 
            scale=20,        
            pivot='middle',
            width=0.001,        # 减小箭头宽度
            linewidth=0.5       # 减小线条宽度
        )
        
        ax.set_title(f"Frame {i:04d}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        
        # 保存当前帧为图片
        frame_filename = os.path.join(save_dir, f"frame_{i:04d}.png")
        plt.savefig(frame_filename, dpi=dpi, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        
        print(f"Saved frame {i:04d} -> {frame_filename}")
    except Exception as e:
        print(f"Error processing frame {i:04d}: {e}")

def main():
    # ---------------------
    #   参数可按需修改
    # ---------------------
    nx = 2048
    ny = 4096
    n_frames = 1024  # 总帧数
    skip = 32        # 矢量箭头采样间隔
    dpi = 100        # 图像分辨率
    
    # 输出帧的保存目录
    save_dir = "frames"
    os.makedirs(save_dir, exist_ok=True)
    
    # 生成所有帧的索引
    frame_indices = list(range(n_frames))
    
    # 使用部分函数传递固定参数
    generate_frame_partial = partial(
        generate_frame, 
        nx=nx, 
        ny=ny, 
        skip=skip, 
        save_dir=save_dir, 
        dpi=dpi
    )
    
    # 确定使用的进程数，通常为 CPU 核心数
    num_processes = cpu_count()
    print(f"Using {num_processes} processes for frame generation.")
    
    # 创建进程池并并行处理帧生成
    with Pool(processes=num_processes) as pool:
        pool.map(generate_frame_partial, frame_indices)
    
    print("All frames have been generated.")

if __name__ == "__main__":
    main()

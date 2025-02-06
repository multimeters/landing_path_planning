import subprocess
import os

# 提供文件列表（至少10个.py文件的名字）
file_list = [
    "1_gen_mask.py", "2_text_to_mat_chunk_processor3.py", "3_cal_angle.py", "4_direction_spectrum_cal.py", "5_fft_cal.py",
    "6_k_cal.py", "7_per_cal2.py", "8_per_to_heave_rao.py", "9_per_to_pitch_rao.py", "10_conv_cal.py", "11_cvar_cal.py", "12_uv_format.py"
]

# 获取当前工作目录
current_directory = os.getcwd()

# 逐个执行文件
for file in file_list:
    file_path = os.path.join(current_directory, file)
    
    if os.path.exists(file_path):
        print(f"正在执行：{file_path}")
        subprocess.run(['python', file_path])  # 使用 subprocess 执行文件
    else:
        print(f"文件 {file} 不存在，跳过。")

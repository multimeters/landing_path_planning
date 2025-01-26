import subprocess
import time

# 设置不同的 exploration 参数值
exploration_values = [0.5, 0.4, 0.3, 0.2, 0.1]
initial_lambda = 300.0
initial_lambda_samples = 2000000

# MATLAB 脚本路径
matlab_script = 'plot_v1.m'

# 循环执行命令
for exploration in exploration_values:


    print(f"执行 MATLAB 脚本: {matlab_script}")
    matlab_command = [
        'matlab', '-nodisplay', '-nosplash', '-r', f"run('{matlab_script}');exit;"
    ]
    
    # 执行 MATLAB 脚本
    result = subprocess.run(matlab_command, capture_output=True, text=True)

    # 输出 MATLAB 脚本执行的结果
    if result.returncode == 0:
        print(f"MATLAB 脚本执行成功，输出：\n{result.stdout}")
    else:
        print(f"MATLAB 脚本执行失败，错误信息：\n{result.stderr}")
    # 构建命令行来执行 ./your_program
    command = [
        '../build/VFRRTPlanner',  # 可执行程序
        '-e', str(exploration),  # exploration 参数
        '-l', str(initial_lambda),  # initial_lambda 参数
        '-s', str(initial_lambda_samples)  # initial_lambda_samples 参数
    ]
    
    # 执行 ./your_program
    print(f"执行: {command}")
    result = subprocess.run(command, capture_output=True, text=True)

    # 输出程序的标准输出和标准错误
    if result.returncode == 0:
        print(f"程序成功执行，输出：\n{result.stdout}")
    else:
        print(f"程序执行失败，错误信息：\n{result.stderr}")

    # 执行 MATLAB 脚本 plot_v1.m
    print(f"执行 MATLAB 脚本: {matlab_script}")
    matlab_command = [
        'matlab', '-nodisplay', '-nosplash', '-r', f"run('{matlab_script}');exit;"
    ]
    
    # 执行 MATLAB 脚本
    result = subprocess.run(matlab_command, capture_output=True, text=True)

    # 输出 MATLAB 脚本执行的结果
    if result.returncode == 0:
        print(f"MATLAB 脚本执行成功，输出：\n{result.stdout}")
    else:
        print(f"MATLAB 脚本执行失败，错误信息：\n{result.stderr}")

    # 可选：根据需要添加延时
    time.sleep(1)  # 延时 1 秒，可以根据需求调整

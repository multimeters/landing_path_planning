import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import numpy as np
import matplotlib.pyplot as plt

# 读取矢量场数据
vector_field_data = np.loadtxt('/home/ganzhi/VFRRT_Planner/vector_Field.txt')
x = vector_field_data[:, 0]
y = vector_field_data[:, 1]
u = vector_field_data[:, 2]
v = vector_field_data[:, 3]

# 读取路径数据
path_data = np.loadtxt('/home/ganzhi/VFRRT_Planner/path_exploration_0.100000_lambda_1.100000_samples_1000000_cost_72.196207.txt')
path_x = path_data[:, 0]
path_y = path_data[:, 1]

# 创建图形和轴
plt.figure(figsize=(10, 8))

# 绘制矢量场的quiver图
quiver = plt.quiver(x, y, u, v, color='gray', alpha=0.7, label='Vector Field')

# 绘制路径
path, = plt.plot(path_x, path_y, 'r', linewidth=2, label='Path')

# 增大起点和终点的标记
start = plt.scatter(path_x[0], path_y[0], s=100, c='g', marker='o', label='Start', zorder=5)
goal = plt.scatter(path_x[-1], path_y[-1], s=100, c='b', marker='o', label='Goal', zorder=5)

# 标注起点和终点
plt.text(path_x[0], path_y[0], ' Start', verticalalignment='bottom',
         horizontalalignment='right', fontsize=12, color='r', fontweight='bold')
plt.text(path_x[-1], path_y[-1], ' Goal', verticalalignment='bottom',
         horizontalalignment='right', fontsize=12, color='b', fontweight='bold')

# 添加标签和标题
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Vector Field and Path')

# 添加图例
plt.legend()

# 添加网格（可选）
plt.grid(True, linestyle='--', alpha=0.5)

# 调整布局
plt.tight_layout()

# 保存图形到文件
plt.savefig('vector_field_and_path.png', dpi=300)

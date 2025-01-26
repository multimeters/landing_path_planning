import matplotlib.pyplot as plt
import numpy as np

# 读取向量场数据
vector_field = np.loadtxt('data/vector_field.txt')
x_vf = vector_field[:,0]
y_vf = vector_field[:,1]
u = vector_field[:,2]
v = vector_field[:,3]

# 读取路径数据
path = np.loadtxt('data/path.txt')
x_path = path[:,0]
y_path = path[:,1]

# 绘制向量场
plt.figure(figsize=(10,10))
plt.quiver(x_vf, y_vf, u, v, color='lightgray', alpha=0.9)

# 绘制路径
plt.plot(x_path, y_path, color='red', linewidth=2, label='Path')

# 绘制起点和终点
plt.plot(x_path[0], y_path[0], 'go', label='Start')
plt.plot(x_path[-1], y_path[-1], 'bo', label='Goal')

# 绘制障碍物
circle = plt.Circle((5,5), 1.0, color='black', fill=True)
plt.gca().add_patch(circle)

plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Vector Field and Planned Path')
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

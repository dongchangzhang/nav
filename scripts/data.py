import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x, y, z = [], [], []

with open('../records/data/1.txt') as f:
    for line in f:

        xyz = line.strip().split(' ')
        x.append(float(xyz[0]))
        y.append(float(xyz[1]))
        z.append(float(xyz[2]))

print(x)
print(y)
print(z)

# 绘制散点图
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(x, y, z)
 
 
# 添加坐标轴(顺序是Z, Y, X)
ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red'})
ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
plt.show()

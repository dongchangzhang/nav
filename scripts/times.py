import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x, y = [], []

with open('../build/times.txt') as f:
    for line in f:

        xy = line.strip().split(' ')
        x.append(float(xy[0]))
        y.append(float(xy[1]))

print(x)
print(y)

plt.plot(x, y)

 
plt.show()

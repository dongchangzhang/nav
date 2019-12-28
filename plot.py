import math
import numpy as np
import matplotlib.pyplot as plt


id = 81
error_file = 'records/errors/%d.txt' % id

X, Y, Z, th, wo, ka, f, x0, y0 = [], [], [], [], [], [], [], [], []

with open(error_file) as file:
    for line in file:
        elements = line.split(' ')
        X.append(float(elements[0]))
        Y.append(float(elements[1]))
        Z.append(float(elements[2]))

        th.append(float(elements[3]) * 180 / math.pi)
        wo.append(float(elements[4]) * 180 / math.pi)
        ka.append(float(elements[5]) * 180 / math.pi)

        f.append(float(elements[6]) * 1e3)
        x0.append(float(elements[7]) * 1e6)
        y0.append(float(elements[8]) * 1e6)

plt.subplot(331)
plt.title('Xs')
plt.xlabel('iterations')
plt.ylabel('error (m)')
plt.plot(range(len(X)), X)

plt.subplot(332)
plt.title('Ys')
plt.xlabel('iterations')
plt.ylabel('error (m)')
plt.plot(range(len(Y)), Y)

plt.subplot(333)
plt.title('Zs')
plt.xlabel('iterations')
plt.ylabel('error (m)')
plt.plot(range(len(Z)), Z)

plt.subplot(334)
plt.title('th')
plt.xlabel('iterations')
plt.ylabel('error: (degrees)')
plt.plot(range(len(th)), th)

plt.subplot(335)
plt.title('wo')
plt.xlabel('iterations')
plt.ylabel('error (degrees)')
plt.plot(range(len(wo)), wo)

plt.subplot(336)
plt.title('ka')
plt.xlabel('iterations')
plt.ylabel('error (degrees)')
plt.plot(range(len(ka)), ka)

plt.subplot(337)
plt.title('f')
plt.xlabel('iterations')
plt.ylabel('error (mm)')
plt.plot(range(len(f)), f)

plt.subplot(338)
plt.title('x0')
plt.xlabel('iterations')
plt.ylabel('error (um)')
plt.plot(range(len(x0)), x0)

plt.subplot(339)
plt.title('y0')
plt.xlabel('iterations')
plt.ylabel('error (um)')
plt.plot(range(len(y0)), y0)

plt.tight_layout()
plt.show()

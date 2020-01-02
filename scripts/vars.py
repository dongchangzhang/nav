import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    exit(-1)

id = int(sys.argv[1])

vars_file = '../records/vars/%d.txt' % id
vars_history_file = '../records/vars_history/%d.txt' % id

X, Y, Z, th, wo, ka, f, x0, y0 = [], [], [], [], [], [], [], [], []

real_vars, base_vars = [], []
with open(vars_file) as file:
    lines = [line.strip().split(' ') for line in list(file)]
    real_vars = [float(x) for x in lines[0]]
    base_vars = [float(x) for x in lines[1]]
    solved = [float(x) for x in lines[2]]
    print('real', real_vars)
    print('base', base_vars)
    print('solved', solved)

with open(vars_history_file) as file:
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
plt.ylabel('Xs-value (m)')
plt.plot(range(len(X)), X)
plt.plot(range(len(X)), [real_vars[0]] * len(X), color='red')
plt.scatter(0, base_vars[0], s=20, c='#FF0000', alpha=0.4)

plt.subplot(332)
plt.title('Ys')
plt.xlabel('iterations')
plt.ylabel('Ys-value (m)')
plt.plot(range(len(Y)), Y)
plt.plot(range(len(Y)), [real_vars[1]] * len(Y), color='red')
plt.scatter(0, base_vars[1], s=20, c='#FF0000', alpha=0.4)

plt.subplot(333)
plt.title('Zs')
plt.xlabel('iterations')
plt.ylabel('Zs-value (m)')
plt.plot(range(len(Z)), Z)
plt.plot(range(len(Z)), [real_vars[2]] * len(Z), color='red')
plt.scatter(0, base_vars[2], s=20, c='#FF0000', alpha=0.4)

plt.subplot(334)
plt.title('th')
plt.xlabel('iterations')
plt.ylabel('th-value (degrees)')
plt.plot(range(len(th)), th)
plt.plot(range(len(th)), [math.degrees(real_vars[3])] * len(th), color='red')
plt.scatter(0, math.degrees(base_vars[3]), s=20, c='#FF0000', alpha=0.4)

plt.subplot(335)
plt.title('wo')
plt.xlabel('iterations')
plt.ylabel('wo-value (degrees)')
plt.plot(range(len(wo)), wo)
plt.plot(range(len(wo)), [math.degrees(real_vars[4])] * len(wo), color='red')
plt.scatter(0, math.degrees(base_vars[4]), s=20, c='#FF0000', alpha=0.4)

plt.subplot(336)
plt.title('ka')
plt.xlabel('iterations')
plt.ylabel('ka-value (degrees)')
plt.plot(range(len(ka)), ka)
plt.plot(range(len(ka)), [math.degrees(real_vars[5])] * len(ka), color='red')
plt.scatter(0, math.degrees(base_vars[5]), s=20, c='#FF0000', alpha=0.4)

plt.subplot(337)
plt.title('f')
plt.xlabel('iterations')
plt.ylabel('f-value (mm)')
plt.plot(range(len(f)), f)
plt.plot(range(len(f)), [real_vars[6] * 1e3] * len(f), color='red')
plt.scatter(0, base_vars[6] * 1e3, s=20, c='#FF0000', alpha=0.4)

plt.subplot(338)
plt.title('x0')
plt.xlabel('iterations')
plt.ylabel('x0-value (um)')
plt.plot(range(len(x0)), x0)
plt.plot(range(len(x0)), [real_vars[7] * 1e6] * len(x0), color='red')
plt.scatter(0, base_vars[7] * 1e6, s=20, c='#FF0000', alpha=0.4)

plt.subplot(339)
plt.title('y0')
plt.xlabel('iterations')
plt.ylabel('y0-value (um)')
plt.plot(range(len(y0)), y0)
plt.plot(range(len(y0)), [real_vars[8] * 1e6] * len(y0), color='red')
plt.scatter(0, base_vars[8] * 1e6, s=20, c='#FF0000', alpha=0.4)

plt.tight_layout()
plt.show()

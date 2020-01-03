import sys
import math
import numpy as np
import matplotlib.pyplot as plt


class Element:
    def __init__(self, real, solved):
        self.real = real
        self.solved = [solved]
        self.len = math.sqrt(real[0] ** 2 + real[1] ** 2 + real[2] ** 2)
    
    def add(self, solved):
        self.solved.append(solved)

class Data:
    def __init__(self):
        self.data = {}
        self.values = []

        self.sorted = False

    def add(self, real, solved):
        key = '%d#%d#%d' % (real[0], real[1], real[2])
        if key in self.data:
            self.data[key].add(solved)
        else:
            self.data[key] = Element(real, solved)

    def split(self):
        self.values = list(self.data.values())
        self.values.sort(key=lambda x: x.len)
        self.sorted = True
    
    def get_value(self, idx):
        if not self.sorted:
            self.split()
        lens = [x.len for x in self.values]
        reals = [x.real[idx] for x in self.values]
        boxes = [[x[idx] for x in y.solved] for y in self.values]
        return lens, reals, boxes


if len(sys.argv) != 3:
    print('need 2 args!')
    exit(-1)

id_b = int(sys.argv[1])
id_e = int(sys.argv[2])

if id_b < 0 or id_e < 0 or id_b >= id_e:
    print('error args!')
    exit(-1)

data = Data()

for idx in range(id_b, id_e + 1):
    # error_file = 'records/errors/%d.txt' % id
    vars_file = '../records/vars/%d.txt' % idx
    real_vars, base_vars, solved = [], [], []
    with open(vars_file) as file:
        lines = [line.strip().split(' ') for line in list(file)]
        real_vars = [float(x) for x in lines[0]]
        base_vars = [float(x) for x in lines[1]]
        solved = [float(x) for x in lines[2]]

        real_vars[3] = math.degrees(real_vars[3])
        real_vars[4] = math.degrees(real_vars[4])
        real_vars[5] = math.degrees(real_vars[5])
        real_vars[6] *= 1e3
        real_vars[7] *= 1e6
        real_vars[8] *= 1e6

        solved[3] = math.degrees(solved[3])
        solved[4] = math.degrees(solved[4])
        solved[5] = math.degrees(solved[5])
        solved[6] *= 1e3
        solved[7] *= 1e6
        solved[8] *= 1e6

        data.add(real_vars, solved)

def draw(idx, where, title, xlabel, ylabel):
    print('--', idx, '---')
    x, y, _y = data.get_value(idx)
    plt.subplot(where)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.gca().invert_xaxis()
    plt.ylabel(ylabel)
    # plt.plot(x, y)
    for i in range(len(_y)):
        for vy in _y[i]:
            if abs(vy - y[i]) < 10:
                plt.scatter(x[i], vy - y[i], s=1, c='red')

# draw(0, 111, 'Xs', 'distance (m)', 'value (m)')
# plt.savefig('../output/statistics.png', dpi=500)
# exit(0)

draw(0, 331, 'Xs', 'distance (m)', 'value (m)')
draw(1, 332, 'Ys', 'distance (m)', 'value (m)')
draw(2, 333, 'Zs', 'distance (m)', 'value (m)')

draw(3, 334, 'th', 'distance (m)', 'value (degrees)')
draw(4, 335, 'wo', 'distance (m)', 'value (degrees)')
draw(5, 336, 'ka', 'distance (m)', 'value (degrees)')

draw(6, 337, 'f', 'distance (m)', 'value (mm)')

draw(7, 338, 'x0', 'distance (m)', 'value (um)')
draw(8, 339, 'y0', 'distance (m)', 'value (um)')

plt.tight_layout()
plt.savefig('../output/statistics.png', dpi=500)
#plt.show()

# Import module
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

dt = 0.05
C = 12
K = 0.9
height = 2
grid = 60

old_H = np.zeros([grid, grid])
H = np.ones([grid, grid])
new_H = np.zeros([grid, grid])
sz = 21

# small peak
z = np.linspace(-1,1,sz)
x = np.ones((sz, sz))
for i in range(sz):
    for j in range(sz):
        x[i][j] = z[i]
y = x.T

TMP_H = height * np.exp(-5 * (x ** 2 + y ** 2))

# half sphere
# z = np.linspace(-2.1, 2.1, sz, dtype=np.float)
# x = np.ones((sz, sz))
# for i in range(sz):
#     for j in range(sz):
#         x[i][j] = z[i]
# y = x.T

# TMP_H = height / 7 * np.sqrt(9 - np.power(x, 2) - np.power(y, 2)) + 1
# print(TMP_H)

H[20:20+sz, 20:20+sz] += np.copy(TMP_H)
old_H = np.copy(H)

x = np.arange(grid)
y = np.arange(grid)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
line = ax.plot_surface(X, Y, H)
ax.view_init(azim=210)
plt.xlabel('x')
plt.ylabel('y')
ax.set_zlim(0, 2)
# plt.ion()
# plt.show()
# plt.close()

def onclick(event):
    global H, old_H, new_H
    bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width*fig.dpi, bbox.height*fig.dpi
    print(event.x)
    print(event.y)
    x = -int(grid * (event.x / width))
    y = int(grid * (event.y / height))
    H[x, y] += 10

# cid = fig.canvas.mpl_connect('button_press_event', onclick)

def data(frame):
    global H, old_H, new_H
    for i in range(grid):
        for j in range(grid):
            if i + 1 >= grid:
                i_add_dx = grid - 1
            else:
                i_add_dx = i + 1
            if i - 1 <= 0:
                i_sub_dx = 0
            else:
                i_sub_dx = i - 1
            if j + 1 >= grid:
                j_add_dy = grid - 1
            else:
                j_add_dy = j + 1
            if j - 1 <= 0:
                j_sub_dy = 0
            else:
                j_sub_dy = j - 1
            # D = H[i_add_dx][j] + H[i_sub_dx][j] + H[i][j_add_dy] + H[i][j_sub_dy]
            # new_H[i][j] = (2 * H[i][j] + pow(C, 2) * (D - 4 * H[i][j]) + old_H[i][j] * (K * dt / 2 - 1) ) / (1 + K * dt / 2)
            P = (H[i_add_dx][j] + H[i_sub_dx][j] + H[i][j_add_dy] + H[i][j_sub_dy] - 4 * H[i][j]) * pow(C, 2)
            # P = (H[i_add_dx][j] + H[i_sub_dx][j] + H[i][j_add_dy] + H[i][j_sub_dy] - 4 * H[i][j] + 0.5 * (H[i_add_dx][j_add_dy] + H[i_sub_dx][j_sub_dy] + H[i_add_dx][j_add_dy] + H[i_sub_dx][j_sub_dy] - 4 * H[i][j] ) ) * pow(C, 2)
            A = -(H[i][j] - old_H[i][j]) / dt + P
            new_H[i][j] = A * pow(dt, 2) + 2 * H[i][j] - old_H[i][j]
    old_H = np.copy(H)
    H = np.copy(new_H)
    ax.clear()
    ax.set_zlim(0, 5)
    line = ax.plot_surface(X, Y, H)
    return line

ani = animation.FuncAnimation(fig, data, fargs=(), interval=30, blit=False)

plt.show()
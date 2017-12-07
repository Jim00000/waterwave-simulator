# -*- encoding:utf-8 -*-
"""
Copyright (C) 2017 the team of Jim00000, ActKz and pityYo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import numpy as np
from mayavi import mlab

dt = 0.04
C = 16
K = 0.1
height = 6
grid = 200

old_H = np.zeros([grid, grid], dtype=np.float64)
H = np.ones([grid, grid], dtype=np.float64)
new_H = np.zeros([grid, grid], dtype=np.float64)
sz = 31

# small peak
z = np.linspace(-1,1,sz)
x = np.ones((sz, sz))
for i in range(sz):
    for j in range(sz):
        x[i][j] = z[i]
y = x.T

TMP_H = height * np.exp(-5 * (x ** 2 + y ** 2))
H[20:20+sz, 20:20+sz] += np.copy(TMP_H)
old_H = np.copy(H)

x = np.arange(grid)
y = np.arange(grid)
X, Y = np.meshgrid(x, y)

def update():
    global H, old_H, new_H

    # Centroid
    for i in range(1, grid - 1):
        for j in range(1, grid - 1):
            P = H[i + 1][j] + H[i - 1][j] + H[i][j + 1] + H[i][j - 1] - 4 * H[i][j]
            new_H[i][j] = ( pow(C * dt, 2) * P * 2 + 4 * H[i][j] - old_H[i][j] * (2 - K * dt) ) / (2 + K * dt)

    # Four edges
    for i in range(1, grid - 1):
        P1 = H[i + 1][0] + H[i - 1][0] + H[i][1] - 3 * H[i][0]
        P2 = H[i + 1][grid - 1] + H[i - 1][grid - 1] + H[i][grid - 2] - 3 * H[i][grid - 1]
        P3 = H[0 + 1][i] + H[0][i + 1] + H[0][i - 1] - 3 * H[0][i]
        P4 = H[grid - 2][i] + H[grid - 1][i + 1] + H[grid - 1][i - 1] - 3 * H[grid - 1][i]
        new_H[i][0] = ( pow(C * dt, 2) * P1 * 2 + 4 * H[i][0] - old_H[i][0] * (2 - K * dt) ) / (2 + K * dt)
        new_H[i][grid - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * H[i][grid - 1] - old_H[i][grid - 1] * (2 - K * dt) ) / (2 + K * dt)
        new_H[0][i] = ( pow(C * dt, 2) * P3 * 2 + 4 * H[0][i] - old_H[0][i] * (2 - K * dt) ) / (2 + K * dt)
        new_H[grid - 1][i] = ( pow(C * dt, 2) * P4 * 2 + 4 * H[grid - 1][i] - old_H[grid - 1][i] * (2 - K * dt) ) / (2 + K * dt)

    # Four corners
    P1 = H[1][0] + H[0][0 + 1] - 2 * H[0][0]
    P2 = H[1][grid - 1] + H[0][grid - 2] - 2 * H[0][grid - 1]
    P3 = H[grid - 2][0] + H[grid - 1][1] - 2 * H[grid - 1][0]
    P4 = H[grid - 2][grid - 1] + H[grid - 1][grid - 2] - 2 * H[grid - 1][grid - 1]
    new_H[0][0] = ( pow(C * dt, 2) * P1 * 2 + 4 * H[0][0] - old_H[0][0] * (2 - K * dt) ) / (2 + K * dt)
    new_H[0][grid-1] = ( pow(C * dt, 2) * P2 * 2 + 4 * H[0][grid-1] - old_H[0][grid-1] * (2 - K * dt) ) / (2 + K * dt)
    new_H[grid-1][0] = ( pow(C * dt, 2) * P3 * 2 + 4 * H[grid-1][0] - old_H[grid-1][0] * (2 - K * dt) ) / (2 + K * dt)
    new_H[grid - 1][grid - 1] = ( pow(C * dt, 2) * P4 * 2 + 4 * H[grid - 1][grid - 1] - old_H[grid - 1][grid - 1] * (2 - K * dt) ) / (2 + K * dt)

    old_H = np.copy(H)
    H = np.copy(new_H)

plt = mlab.surf(H, warp_scale='auto', colormap=u'ocean')

@mlab.animate(delay=10)
def animation():
    f = mlab.gcf()
    while True:
        update()
        plt.mlab_source.set(scalars=H)
        f.scene.render()
        yield


animation()
mlab.title('sequential in Python')
mlab.show()
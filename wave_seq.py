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
from wave_equation import sequential_update

dt = 0.04
C = 16
K = 0
height = 6
grid = 100

old_H = np.zeros([grid, grid], dtype=np.float64)
H = np.ones([grid, grid], dtype=np.float64)
new_H = np.zeros([grid, grid], dtype=np.float64)
sz = 5

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

def update():
    global H, old_H, new_H
    H, old_H, new_H = sequential_update(H, old_H, new_H, grid, grid, C, K, dt) 

plt = mlab.surf(H, warp_scale='auto', colormap=u'ocean')

@mlab.animate(delay=10)
def animation():
    f = mlab.gcf()
    while True:
        update()
        plt.mlab_source.set(scalars=H)
        f.scene.render()
        yield
raw_input()
animation()
# mlab.title('sequential in C')
mlab.show()
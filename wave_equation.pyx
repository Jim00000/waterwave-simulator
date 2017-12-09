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
cimport numpy as np

cdef extern from "wave_equation_api.h":
    void c_sequential_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)
    void c_openmp_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)
    void c_threadpool_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)

cdef extern from "wave_equation_cuda_api.h":
    void c_cuda_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)

def sequential_update(np.ndarray[double, ndim=2, mode="c"] data not None, np.ndarray[double, ndim=2, mode="c"] olddata not None, np.ndarray[double, ndim=2, mode="c"] newdata not None, row_size, col_size, C, K, dt):
    c_sequential_update(&data[0, 0], &olddata[0, 0], &newdata[0, 0], row_size, col_size, C, K, dt)
    return newdata, data, olddata

def openmp_update(np.ndarray[double, ndim=2, mode="c"] data not None, np.ndarray[double, ndim=2, mode="c"] olddata not None, np.ndarray[double, ndim=2, mode="c"] newdata not None, row_size, col_size, C, K, dt):
    c_openmp_update(&data[0, 0], &olddata[0, 0], &newdata[0, 0], row_size, col_size, C, K, dt)
    return newdata, data, olddata

def threadpool_update(np.ndarray[double, ndim=2, mode="c"] data not None, np.ndarray[double, ndim=2, mode="c"] olddata not None, np.ndarray[double, ndim=2, mode="c"] newdata not None, row_size, col_size, C, K, dt):
    c_threadpool_update(&data[0, 0], &olddata[0, 0], &newdata[0, 0], row_size, col_size, C, K, dt)
    return newdata, data, olddata

def cuda_update(np.ndarray[double, ndim=2, mode="c"] data not None, np.ndarray[double, ndim=2, mode="c"] olddata not None, np.ndarray[double, ndim=2, mode="c"] newdata not None, row_size, col_size, C, K, dt):
    c_cuda_update(&data[0, 0], &olddata[0, 0], &newdata[0, 0], row_size, col_size, C, K, dt)
    return newdata, data, olddata

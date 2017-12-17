/*
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
*/

/**
 * @file        wave_equation_cuda_api.cu
 * @author      Jim00000
 * @date        12.8.2017
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "wave_equation_cuda_api.h"

__global__ void _cuda_update_kernal_(double* d_data, double* d_olddata, double* d_newdata, int row_size, int col_size, double C, double K, double dt);

__global__ void _cuda_update_kernal_(double* d_data, double* d_olddata, double* d_newdata, int row_size, int col_size, double C, double K, double dt)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int idy = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = idy, j = idx;

    double potential = d_data[(i + 1) * col_size + j] + d_data[(i - 1) * col_size + j] + d_data[i * col_size + j + 1] +
    d_data[i * col_size + j - 1] - 4 * d_data[i * col_size + j];
    
    d_newdata[i * col_size + j] = ( pow(C * dt, 2) * potential * 2 + 4 * d_data[i * col_size + j] - d_olddata[i * col_size + j] *
    (2 - K * dt) ) / (2 + K * dt);
}

void c_cuda_update(double* data, double* olddata, double* newdata, int row_size, int col_size, double C, double K, double dt)
{
    const static int ARRAY_SIZE = row_size * col_size;
    const static int ARRAY_BYTES = ARRAY_SIZE * sizeof(double);     
    static bool initialized = false;

    // Declare GPU memory pointers
    static double* d_data;
    static double* d_olddata;
    static double* d_newdata;

    if(initialized == false) {
        initialized = true;
        // Allocate GPU memory
        cudaMalloc((void**) &d_data, ARRAY_BYTES);
        cudaMalloc((void**) &d_olddata, ARRAY_BYTES);
        cudaMalloc((void**) &d_newdata, ARRAY_BYTES);
    }


    // Transfer memory to GPU memory
    cudaMemcpy(d_data, data, ARRAY_BYTES, cudaMemcpyHostToDevice);
    cudaMemcpy(d_olddata, olddata, ARRAY_BYTES, cudaMemcpyHostToDevice);

    const static int THREADS_COUNT = 16;
    dim3 threads(THREADS_COUNT, THREADS_COUNT);
    dim3 blocks(col_size / threads.x + 1, row_size / threads.y + 1);

    // Launch the kernel
    _cuda_update_kernal_<<<blocks, threads>>>(d_data, d_olddata, d_newdata, row_size, col_size, C, K, dt);
    cudaDeviceSynchronize();

    // Transfer GPU memory to host memory
    cudaMemcpy(data, d_data, ARRAY_BYTES, cudaMemcpyDeviceToHost);
    cudaMemcpy(newdata, d_newdata, ARRAY_BYTES, cudaMemcpyDeviceToHost);

    // Free GPU memory
    // cudaFree(d_data);
    // cudaFree(d_olddata);
    // cudaFree(d_newdata);

    // Four edges
    #pragma omp parallel for shared(row_size, col_size, newdata, data, olddata, C, K, dt)
    for(int i = 1; i < row_size - 1; i++) {
        double P1 = data[(i + 1) * col_size] + data[(i - 1) * col_size] + data[i * col_size + 1] - 3 * data[i * col_size];
        double P2 = data[(i + 1) * col_size + col_size - 1] + data[(i - 1) * col_size + col_size - 1] +
                    data[i * col_size + col_size - 2] - 3 * data[i * col_size + col_size - 1];
        double P3 = data[col_size + i] + data[i + 1] + data[i - 1] - 3 * data[i];
        double P4 = data[(row_size - 2) * col_size + i] + data[(row_size - 1) * col_size + i + 1] +
                    data[(row_size - 1) * col_size + i - 1] - 3 * data[(row_size - 1) * col_size + i];
        newdata[i * col_size] = ( pow(C * dt, 2) * P1 * 2 + 4 * data[i * col_size] - olddata[i * col_size] *
                                    (2 - K * dt) ) / (2 + K * dt);
        newdata[i * col_size + col_size - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * data[i * col_size + col_size - 1] -
                                                olddata[i * col_size + col_size - 1] * (2 - K * dt) ) / (2 + K * dt);
        newdata[i] = ( pow(C * dt, 2) * P3 * 2 + 4 * data[i] - olddata[i] * (2 - K * dt) ) / (2 + K * dt);
        newdata[(row_size - 1) * col_size + i] = ( pow(C * dt, 2) * P4 * 2 + 4 * data[(row_size - 1) * col_size + i] -
                olddata[(row_size - 1) * col_size + i] * (2 - K * dt) ) / (2 + K * dt);
    }

    // Four corners
    double P1 = data[col_size] + data[1] - 2 * data[0];
    double P2 = data[col_size + col_size - 1] + data[col_size - 2] - 2 * data[col_size - 1];
    double P3 = data[(row_size - 2) * col_size] + data[(row_size - 1) * col_size + 1] - 2 * data[(row_size - 1) * col_size];
    double P4 = data[(row_size - 2) * col_size + col_size - 1] + data[(row_size - 1) * col_size + col_size - 2] - 2 *
                data[(row_size - 1) * col_size + col_size - 1];
    newdata[0] = ( pow(C * dt, 2) * P1 * 2 + 4 * data[0] - olddata[0] * (2 - K * dt) ) / (2 + K * dt);
    newdata[col_size - 1] = ( pow(C * dt, 2) * P2 * 2 + 4 * data[col_size - 1] - olddata[col_size - 1] * (2 - K * dt) ) /
                            (2 + K * dt);
    newdata[(row_size - 1) * col_size] = ( pow(C * dt, 2) * P3 * 2 + 4 * data[(row_size - 1) * col_size] - olddata[(row_size - 1)
                                            * col_size] * (2 - K * dt) ) / (2 + K * dt);
    newdata[(row_size - 1) * col_size + col_size - 1] = ( pow(C * dt, 2) * P4 * 2 +
            4 * data[(row_size - 1) * col_size + col_size - 1] - olddata[(row_size - 1) * col_size + col_size - 1] * (2 - K * dt) )
            / (2 + K * dt);

}

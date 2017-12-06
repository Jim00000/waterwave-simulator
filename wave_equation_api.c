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
 * @file        wave_equation_api.c
 * @author      Jim00000
 * @date        12.6.2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "wave_equation_api.h"

#define C 12
#define dt 0.05

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

inline int select_add_index(const int index, const int max_edge) __attribute__((always_inline));
inline int select_sub_index(const int index) __attribute__((always_inline));

void c_sequential_update(double* data, double* olddata, double* newdata, int row_size, int col_size)
{
    // Check that data is not a null pointer
    if(data == NULL) {
        fprintf(stderr, "[ERROR] 2D array 'data' is null \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that row_size is a positive number
    if(row_size <= 0) {
        fprintf(stderr, "[ERROR] row_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that col_size is a positive number
    if(col_size <= 0) {
        fprintf(stderr, "[ERROR] col_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < row_size; i++) {
        for(int j = 0; j < col_size; j++) {
            int i_add_dx = select_add_index(i, row_size);
            int j_add_dy = select_add_index(j, col_size);
            int i_sub_dx = select_sub_index(i);
            int j_sub_dy = select_sub_index(j);
            double potential = (data[i_add_dx * col_size + j] + data[i_sub_dx * col_size + j] +
                                data[i * col_size + j_add_dy] + data[i * col_size + j_sub_dy] -
                                4 * data[i * col_size + j]) * pow(C, 2);
            double acceleration = -(data[i * col_size + j] - olddata[i * col_size + j]) / dt + potential;
            double newvalue = acceleration * pow(dt, 2) + 2 * data[i * col_size + j] - olddata[i * col_size + j];
            newdata[i * col_size + j] = newvalue;
        }
    }

}

void c_openmp_update(double* data, double* olddata, double* newdata, int row_size, int col_size)
{
    // Check that data is not a null pointer
    if(data == NULL) {
        fprintf(stderr, "[ERROR] 2D array 'data' is null \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that row_size is a positive number
    if(row_size <= 0) {
        fprintf(stderr, "[ERROR] row_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    // Make sure that col_size is a positive number
    if(col_size <= 0) {
        fprintf(stderr, "[ERROR] col_size must be a positive number \n");
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for shared(row_size, col_size, data, olddata, newdata) collapse(2)
    for(int i = 0; i < row_size; i++) {
        for(int j = 0; j < col_size; j++) {
            int i_add_dx = select_add_index(i, row_size);
            int j_add_dy = select_add_index(j, col_size);
            int i_sub_dx = select_sub_index(i);
            int j_sub_dy = select_sub_index(j);
            double potential = (data[i_add_dx * col_size + j] + data[i_sub_dx * col_size + j] +
                                data[i * col_size + j_add_dy] + data[i * col_size + j_sub_dy] -
                                4 * data[i * col_size + j]) * pow(C, 2);
            double acceleration = -(data[i * col_size + j] - olddata[i * col_size + j]) / dt + potential;
            double newvalue = acceleration * pow(dt, 2) + 2 * data[i * col_size + j] - olddata[i * col_size + j];
            newdata[i * col_size + j] = newvalue;
        }
    }

}

int select_add_index(const int index, const int max_edge)
{
    return unlikely(index + 1 >= max_edge) ? max_edge - 1 : index + 1;
}

int select_sub_index(const int index)
{
    return unlikely(index - 1 <= 0) ? 0 : index - 1;
}